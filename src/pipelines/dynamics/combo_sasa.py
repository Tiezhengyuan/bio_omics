'''
requirments:
    - table combo2, combo3, combo4, combo5
example: 
    - parallel -j8 bomics combo_sasa ::: {0..241}
    - bomics combo_sasa 0
functions:
    - calculate SASA and analyze PPI
    - put meta data to table combo_sasa
'''
from src.helper import *
from bioomics import SasaProcess, ProcessPickle, TABLE_COLUMNS

@header
def retrieve_query(params, meta):
    combo = params['combo']
    combo_table = f"combo{combo}"
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        SELECT DISTINCT S.pdb_id, S.combo_id, S.chain_combo, S.combo_pdb, S.chain_pdbs
        FROM label_complex      C
        LEFT JOIN {combo_table} S ON C.combo_id = S.combo_id
        LEFT JOIN meta_pdb      P ON S.pdb_id = P.pdb_id
        WHERE C.combo={combo}
            AND P.pdb_group = {pdb_group}
            AND S.combo_pdb IS NOT NULL
            AND NOT EXISTS (
                SELECT 1 FROM {table_name} x
                WHERE x.combo_id = C.combo_id
            )
    """
    rows = QueryComplex(params['verbose']).list_data(query)
    meta['rows'] = len(rows)
    print(f"Group {pdb_group}: Collect {combo_table} chains from {len(rows)} PDB ...")
    for row in rows:
        row['combo_pdb_file'] = params['output_pdb_dir'] / row['combo_pdb']
        chain_pdbs = row['chain_pdbs'].split(',')
        row['chain_pdb_files'] = [params['output_pdb_dir'] / f for f in chain_pdbs]
        yield row

def cal_sasa(row_iter, params, meta):
    # pdb_id, combo_id, combo_pdb_file, chain_pdb_files
    for row in row_iter:
        pdb_id = row['pdb_id']
        # output file
        endpoint = Path(params['table_name']) / pdb_id[:2] / pdb_id
        outdir = params['simulate_dir'] / endpoint
        outdir.mkdir(parents=True, exist_ok=True)
        file_name =  f"{row['combo_id']}.pkl"
        pfile = outdir / file_name

        if params['overwrite'] or not pfile.is_file():
            try:
                pdb_files_pool = [
                    (row['combo_pdb_file'], row['chain_pdb_files']),
                ]
                SasaProcess(verbose=params['verbose']).process(pdb_files_pool, pfile)
                meta['sasa_succeed'] += 1
            except Exception as e:
                print(e)
                meta['sasa_fail'] += 1
        else:
            meta['skip'] += 1
        
        # build record
        item = {
            'pdb_id': pdb_id,
            'combo_id': row['combo_id'],
        }
        if pfile.is_file():
            item['relative_pkl'] = str(endpoint / file_name)
            data = ProcessPickle(pfile, params['verbose']).load_pickle([])
            if data:
                item['delta_total_sasa'] = data[0].get('delta_total_sasa')
        rec = [item.get(col) for col in TABLE_COLUMNS[params['table_name']]]
        yield tuple(rec)
        break

if __name__ == "__main__":
    params.update({
        'overwrite': False,
        'chunk_size': 20,
        'table_name': 'combo_sasa',
    })
    # clean null data in the table
    validate_data('relative_pkl', params, meta, params['simulate_dir'])

    
    for combo in (2,3,4,5):
        params['combo'] = combo
        
        # calculate delta total SASA
        rows = retrieve_query(params, meta)
        record_iter = cal_sasa(rows, params, meta)

        bc = BuildComplex(params['verbose'], params['chunk_size'])
        bc.insert_batch_records(record_iter, params['table_name'])
        break

    footer(meta)
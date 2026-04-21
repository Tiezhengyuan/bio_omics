'''
requirments:
    - table combo2
example: 
    - parallel -j8 python app.py combo2_sasa ::: {0..237}
    - python app.py combo2_sasa 0
functions:
    - calculate SASA and analyze PPI
    - put meta data to table combo2_sasa
'''
from src.helper import *
from bioomics import SasaProcess

@header
def retrieve_query(params, meta):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        WITH temp AS (
            SELECT pdb_id FROM pdb_complex WHERE combo=2
            EXCEPT
            SELECT pdb_id FROM {table_name}
        )
        SELECT S.pdb_id, S.combo_id, S.chain_combo, S.combo_pdb, S.chain_pdbs
        FROM combo2 S
        LEFT JOIN meta_pdb P ON S.pdb_id = P.pdb_id
        WHERE P.pdb_group = {pdb_group}
            AND S.combo_pdb IS NOT NULL
            AND S.pdb_id IN (
                SELECT pdb_id FROM temp
            )
    ;"""
    rows = QueryComplex(params['verbose']).list_data(query)
    meta['rows'] = len(rows)
    print(f"Group {pdb_group}: Collect combo2 chains from {len(rows)} PDB ...")
    data_dir = params['data_dir']
    for row in rows:
        row['combo_pdb_file'] = os.path.join(data_dir, row['combo_pdb'])
        row['chain_pdb_files'] = [os.path.join(data_dir, f) for f in row['chain_pdbs'].split(',')]
        yield row

def cal_sasa(rows, params, meta):
    # pdb_id, combo_id, combo_pdb_file, chain_pdb_files
    for row in rows:
        pdb_id = row['pdb_id']
        endpoint = os.path.join(params['table_name'], pdb_id[:2], pdb_id)
        outdir = os.path.join(params['data_dir'], endpoint)
        Path(outdir).mkdir(parents=True, exist_ok=True)

        file_name =  f"{row['combo_id']}.pkl"
        pfile = os.path.join(outdir, file_name)
        if params['overwrite'] or (not os.path.isfile(pfile)):
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
        delta_sasa, relative_path = None, None
        if os.path.isfile(pfile):
            relative_path = os.path.join(endpoint, file_name)
            data = ProcessPickle(pfile, params['verbose']).load_pickle({})
            delta_sasa = data[0].get('delta_total_sasa')
        yield (pdb_id, row['combo_id'], delta_sasa, relative_path)

if __name__ == "__main__":
    params.update({
        'overwrite': False,
        'chunk_size': 20,
        'table_name': 'combo2_sasa',
    })
    # calculate delta total SASA
    rows = retrieve_query(params, meta)
    record_iter = cal_sasa(rows, params, meta)

    bc = BuildComplex(record_iter, params['verbose'], params['chunk_size'])
    bc.insert_batch_records(params['table_name'])

    validate_data('relative_pkl', params, meta)
    footer(meta)

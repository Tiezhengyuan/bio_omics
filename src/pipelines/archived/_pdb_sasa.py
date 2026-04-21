'''
requirments:
    - PDB raw data
    - table meta_pdb
example: 
    - parallel -j8 python src/pipelines/pdb_sasa.py ::: {0..117}
    - python src/pipelines/pdb_sasa.py 0
functions:
    - calculate SASA and analyze PPI
    - put meta data to table pdb_sasa
'''
from helper import *
from bioomics import SasaProcess

@header
def retrieve_query(params, meta):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        SELECT S.pdb_id, S.combo_id, S.combo_pdb, S.chain_pdbs
        FROM meta_combo2   S
        LEFT JOIN meta_pdb P ON S.pdb_id = P.pdb_id
        WHERE P.pdb_group = {pdb_group}
            AND P.relative_path IS NOT NULL
            AND S.combo_id IS NOT NULL
            AND S.combo_pdb IS NOT NULL
            AND S.chain_pdbs IS NOT NULL
            AND S.combo_id NOT IN(
                SELECT combo_id FROM {table_name}
            )
            AND S.pdb_id IN (
                SELECT pdb_id from view_antibody
            )
    ;"""
    rows = QueryComplex(params['verbose']).list_data(query)
    meta['rows'] = len(rows)
    print(f"Group {pdb_group}: Collect combo2 chains from {len(rows)} PDB ...")
    for row in rows:
        combo_pdb_file = os.path.join(params['data_dir'], row['combo_pdb'])
        if os.path.isfile(combo_pdb_file):
            row['combo_pdb_file'] = combo_pdb_file
            chain_pdbs = row['chain_pdbs'].split(',')
            if chain_pdbs:
                chain_pdb_files = []
                for relative_path in chain_pdbs:
                    path = os.path.join(params['data_dir'], relative_path)
                    if os.path.isfile(path):
                        chain_pdb_files.append(path)
                if len(chain_pdb_files) == len(chain_pdbs):
                    row['chain_pdb_files'] = chain_pdb_files
                    yield row
            else:
                meta['miss_chain_pdb'] += 1
        else:
            meta['no_combo_pdb'] += 1

def cal_sasa(rows, params, meta):
    for row in rows:
        pdb_id = row['pdb_id']
        endpoint = os.path.join(params['table_name'], pdb_id[:2], pdb_id)
        outdir = os.path.join(params['data_dir'], endpoint)
        Dir(outdir).init_dir()

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
        relative_path = os.path.join(endpoint, file_name) if os.path.isfile(pfile) else None
        yield (pdb_id, row['combo_id'], relative_path)

if __name__ == "__main__":
    params = {
        'pdb_group': sys.argv[1],
        'verbose': True,
        'overwrite': False,
        'chunk_size': 20,
        'data_dir': os.getenv('data_dir'),
        'table_name': 'pdb_sasa',
    }
    # calculate delta total SASA
    rows = retrieve_query(params, meta)
    record_iter = cal_sasa(rows, params, meta)

    bc = BuildComplex(record_iter, params['verbose'], params['chunk_size'])
    bc.insert_batch_records(params['table_name'])

    validate_data('relative_pkl', params, meta)
    footer(meta)

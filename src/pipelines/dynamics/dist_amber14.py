'''
requirments:
    - PDB raw data
    - table relaxed_amber14_combo2, only calculate distance of two chains 
example: 
    - parallel -j8 python src/pipelines/dist_amber14.py ::: {0..117}
    - python src/pipelines/dist_amber14.py 0
functions:
    - calculate distance and analyze PPI
    - put meta data to table dist_amber14
'''
from helper import *
from bioomics import DistanceCalculate, ProcessPickle

@header
def retrieve_query(params, meta):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        SELECT S.pdb_id, S.combo_id, S.chain_combo, S.endpoint_prefix
        FROM relaxed_amber14_combo2  S
        LEFT JOIN meta_pdb           P ON S.pdb_id = P.pdb_id
        WHERE P.pdb_group = {pdb_group}
            AND S.combo_id IS NOT NULL
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
        input_pdb = os.path.join(params['data_dir'], row['endpoint_prefix']+'.pdb')
        if os.path.isfile(input_pdb):
            row['input_pdb'] = input_pdb
            yield row

def cal_dist(rows, params, meta):
    for row in rows:
        pdb_id = row['pdb_id']
        endpoint = os.path.join(params['table_name'], pdb_id[:2], pdb_id)
        outdir = os.path.join(params['data_dir'], endpoint)
        Dir(outdir).init_dir()

        file_name =  f"{row['combo_id']}.pkl"
        pfile = os.path.join(outdir, file_name)
        if params['overwrite'] or (not os.path.isfile(pfile)):
            result = DistanceCalculate(row['input_pdb'], params).process(row['chain_combo'])
            if result:
                result.update(row)
                ProcessPickle(pfile, params['verbose']).save_pickle(result)
                meta['dist_succeed'] += 1
            else:
                meta['dist_fail'] += 1
        else:
            meta['skip_cal_dist'] += 1
        
        # build record
        relative_path = os.path.join(endpoint, file_name) if os.path.isfile(pfile) else None
        yield pdb_id, row['combo_id'], relative_path

if __name__ == "__main__":
    params = {
        'pdb_group': sys.argv[1],
        'verbose': True,
        'overwrite': False,
        'chunk_size': 20,
        'data_dir': os.getenv("data_dir"),
        'table_name': 'dist_amber14',
        'dist_cutoff': 10,
        'span': 10,
    }
    pdb_iter = retrieve_query(params, meta)
    # calculate atom distance
    record_iter = cal_dist(pdb_iter, params, meta)

    # put meta data to table
    bc = BuildComplex(record_iter, params['verbose'], params['chunk_size'])
    bc.insert_batch_records(params['table_name'])

    validate_data('relative_pkl', params, meta)
    footer(meta)

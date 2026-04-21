'''
requirments:
    - PDB raw data
    - table combo2, only calculate distance of two chains 
example: 
    - parallel -j8 python app.py combo2_distance ::: {0..237}
    - python app.py combo2_distance 0
functions:
    - calculate distance and analyze PPI
    - put meta data to table combo2_distance
'''
from src.helper import *
from bioomics import DistanceCalculate, DistanceAnalyze

def retrieve_query(params, meta):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        WITH temp AS (
            SELECT pdb_id FROM pdb_complex WHERE combo=2
            EXCEPT
            SELECT pdb_id FROM {table_name}
        )
        SELECT S.pdb_id, S.combo_id, S.chain_combo, S.combo_pdb
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
    return rows

def cal_dist(rows, params, meta) -> tuple:
    for row in rows:
        pdb_id = row['pdb_id']
        endpoint = os.path.join(params['table_name'], pdb_id[:2], pdb_id)
        outdir = os.path.join(params['data_dir'], endpoint)
        Path(outdir).mkdir(parents=True, exist_ok=True)

        file_name =  f"{row['combo_id']}.pkl"
        pfile = os.path.join(outdir, file_name)
        if params['overwrite'] or (not os.path.isfile(pfile)):
            input_pdb = os.path.join(params['data_dir'], row['combo_pdb'])
            result = DistanceCalculate(input_pdb, params).process(row['chain_combo'])
            if result:
                result.update(row)
                ProcessPickle(pfile, params['verbose']).save_pickle(result)
                meta['dist_succeed'] += 1
            else:
                meta['dist_fail'] += 1
        else:
            meta['skip_cal_dist'] += 1
        
        # build record
        rec0 = (pdb_id, row['combo_id'])
        rec1 = (None,) * 4
        if os.path.isfile(pfile):
            relative_path = os.path.join(endpoint, file_name)
            p = DistanceAnalyze({'pickle_file': pfile})
            dist_value, res1, res2 = p.get_distance(False, 0)
            rec1 = (relative_path, dist_value, res1, res2)
        yield rec0 + rec1


if __name__ == "__main__":
    params.update({
        'overwrite': False,
        'chunk_size': 20,
        'table_name': 'combo2_distance',
        'dist_cutoff': 10,
        'span': 10,
    })
    pdb_iter = retrieve_query(params, meta)

    # calculate atom distance
    record_iter = cal_dist(pdb_iter, params, meta)

    # put meta data to table
    bc = BuildComplex(record_iter, params['verbose'], params['chunk_size'])
    bc.insert_batch_records(params['table_name'])

    validate_data('relative_pkl', params, meta)
    footer(meta)

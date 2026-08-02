'''
requirments:
    - PDB raw data
    - table combo2, only calculate distance of two chains 
example: 
    - parallel -j8 bomics combo2_distance ::: {0..241}
    - bomics combo2_distance 0
functions:
    - calculate distance and analyze PPI
    - put meta data to table combo2_distance
'''
from src.helper import *
from bioomics import TABLE_COLUMNS, ProcessPickle, DistanceCalculate, DistanceAnalyze

def retrieve_query(params, meta):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        WITH temp AS (
            SELECT pdb_id FROM label_complex WHERE combo=2
            EXCEPT
            SELECT pdb_id FROM {table_name}
        )
        SELECT S.pdb_id, S.combo_id, S.chain_combo, S.combo_pdb
        FROM combo2        S
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

def cal_dist(rows, params, meta):
    for row in rows:
        pdb_id = row['pdb_id']
        # define output
        endpoint = Path(params['table_name']) / pdb_id[:2] / pdb_id
        outdir = params['simulate_dir'] / endpoint
        outdir.mkdir(parents=True, exist_ok=True)
        file_name =  f"{row['combo_id']}.pkl"
        pfile = outdir / file_name

        if params['overwrite'] or not pfile.is_file():
            input_pdb = params['output_pdb_dir'] / row['combo_pdb']
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
        item = {
            'pdb_id': pdb_id,
            'combo_id': row['combo_id'],
        }
        if pfile.is_file():
            item['relative_pkl'] = str(endpoint / file_name)
            data = ProcessPickle(pfile).load_pickle([])
            dist_value, res1, res2 = DistanceAnalyze.get_distance(data, False, 0)
            item.update({
                'min_dist': dist_value,
                'min_res1': res1,
                'min_res2': res2,
            })
        cols = TABLE_COLUMNS[params['table_name']]
        rec = [item.get(col) for col in cols]
        yield tuple(rec)


if __name__ == "__main__":
    params.update({
        'overwrite': False,
        'chunk_size': 20,
        'table_name': 'combo2_distance',
        'dist_cutoff': 10,
        'span': 10,
    })
    validate_data('relative_pkl', params, meta, params['simulate_dir'])

    # calculate atom distance
    pdb_iter = retrieve_query(params, meta)
    record_iter = cal_dist(pdb_iter, params, meta)

    # put meta data to table
    bc = BuildComplex(params['verbose'], params['chunk_size'])
    bc.insert_batch_records(record_iter, params['table_name'])

    footer(meta)
'''
example: 
    - parallel -j8 python app.py predict_dist ::: {0..10}
    - python app.py predict_dist 0
functions:
    - calculate distance and analyze PPI
    - put meta data to table predict_dist
'''
from src.helper import *
from bioomics import TABLE_COLUMNS, ProcessPickle, DistanceCalculate, DistanceAnalyze

def retrieve(params, meta):
    min_id = int(params['pdb_group']) * 1_000
    max_id = min_id + 1_000
    table_name = params['table_name']
    col_name = params['col_name']
    query = f"""
        SELECT summary_id, predict_id, {col_name} AS relative_pdb
        FROM predict_summary
        WHERE predict_id >= {min_id} AND predict_id < {max_id}
            AND {col_name} IS NOT NULL
            AND summary_id NOT IN (
                SELECT summary_id from {table_name}
            )
    """
    rows = QueryComplex().list_data(query)
    meta['rows'] = len(rows)
    for row in rows:
        pdb_file = os.path.join(params['design_dir'], row['relative_pdb'])
        if os.path.isfile(pdb_file):
            row['pdb_file'] = pdb_file
            yield row

def cal_dist(rows, params, meta):
    cols = TABLE_COLUMNS[params['table_name']]
    for row in rows:
        predict_id = str(row['predict_id'])
        endpoint = os.path.join(params['table_name'], predict_id[:2], predict_id)
        output_dir = os.path.join(params['design_dir'], endpoint)
        Path(output_dir).mkdir(parents=True, exist_ok=True)

        # calculate distance
        file_name = re.sub('pdb$', 'pkl', os.path.basename(row['relative_pdb']))
        pfile = os.path.join(output_dir, file_name)
        if params['overwrite'] or (not os.path.isfile(pfile)):
            res = DistanceCalculate(row['pdb_file'], params)()
            if res:
                ProcessPickle(pfile, params['verbose']).save_pickle(res)
                meta['dist_succeed'] += 1
            else:
                meta['dist_fail'] += 1
        else:
            meta['skip_cal_dist'] += 1
        
        # build record
        if os.path.isfile(pfile):
            data_pool = ProcessPickle(pfile).load_pickle([])
            for data in data_pool:
                dist_value, res1, res2 = DistanceAnalyze.get_distance(data, False, 0)
                item = {
                    'summary_id': row['summary_id'],
                    'relative_pdb': row['relative_pdb'],
                    'relative_pkl': os.path.join(endpoint, file_name),
                    'chain_combo': data.get('chain_combo'),
                    'min_dist': float(dist_value),
                    'min_res1': res1,
                    'min_res2': res2,
                }
                rec = [item.get(col) for col in cols]
                yield tuple(rec)
        else:
            rec = [row.get(col) for col in cols]
            yield tuple(rec)


if __name__ == "__main__":
    params.update({
        'overwrite': False,
        'chunk_size': 20,
        'dist_cutoff': 10,
        'span': 10,
    })

    for name in ('unrelaxed', 'relaxed'):
        params['table_name'] = f'{name}_dist'
        params['col_name'] = f"{name}_pdb"
        pdb_iter = retrieve(params, meta)
        record_iter = cal_dist(pdb_iter, params, meta)

        bc = BuildComplex(params['verbose'], params['chunk_size'])
        bc.insert_batch_records(record_iter, params['table_name'])

    footer(meta)


'''
requirements:
    - table predict
example: 
    - parallel -j8 python app.py predict_summary ::: {0..10}
    - python app.py predict_summary 0
functions:
    - retrieve summary of alphafold prediction to table predict_summary
'''

from src.helper import *
from bioomics import ProcessJson

def parse_data(row, params):
    res = {}
    indir = os.path.join(params['design_dir'], row['endpoint'])
    for file_name in os.listdir(indir):
        path = os.path.join(indir, file_name)
        if 'rank_' in file_name:
            ranking = int(re.findall(r'rank_(\d+)_', file_name)[0])
            if ranking not in res:
                res[ranking] = {'ranking': ranking,}
            relative_path = os.path.join(row['endpoint'], file_name)
            if '_scores_' in file_name and file_name.endswith('.json'):
                data = ProcessJson.load_json(path)
                res[ranking]['avg_plddt'] = np.mean(data['plddt'])
                res[ranking]['avg_ptm'] = np.mean(data['ptm'])
                res[ranking]['max_pae'] = data['max_pae']
                res[ranking]['relative_json'] = relative_path
            elif '_relaxed_' in file_name and file_name.endswith('.pdb'):
                res[ranking]['relaxed_pdb'] = relative_path
            elif '_unrelaxed_' in file_name and file_name.endswith('.pdb'):
                res[ranking]['unrelaxed_pdb'] = relative_path
    return res

@header
def retrieve_rows(params, meta):
    min_id = int(params['pdb_group']) * 1_000
    max_id = min_id + 1_000
    table_name = params['table_name']
    query = f"""
        SELECT relaxed_pdb FROM {table_name}
        WHERE predict_id >= {min_id} AND predict_id < {max_id}
    ;"""
    curr_rows = QueryComplex(params['verbose']).list_data(query)
    pdb_files = [i['relaxed_pdb'] for i in curr_rows]

    query = f"""
        SELECT predict_id, endpoint FROM predict
        WHERE predict_id >= {min_id} AND predict_id < {max_id}
            AND endpoint IS NOT NULL
    ;"""
    rows = QueryComplex(params['verbose']).list_data(query)
    meta['predict'] = len(rows)

    for row in rows:
        ranked = parse_data(row, params)
        for ranking, log in ranked.items():
            relaxed_pdb = log.get('relaxed_pdb')
            if relaxed_pdb and relaxed_pdb not in pdb_files:
                yield (
                    f"{row['predict_id']}_{ranking}",
                    row['predict_id'],
                    ranking,
                    log.get('avg_plddt'),
                    log.get('avg_ptm'),
                    log.get('max_pae'),
                    log.get('relative_json'),
                    relaxed_pdb, 
                    log.get('unrelaxed_pdb'),
                )
        meta['records'] += len(ranked)
 
if __name__ == "__main__":
    params.update({
        'table_name': 'predict_summary',
    })

    record_iter = retrieve_rows(params, meta)
    bc = BuildComplex(params['verbose'])
    bc.insert_batch_records(record_iter, params['table_name'])

    footer(meta)

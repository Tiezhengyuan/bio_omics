'''
requirements:
    - table chain_alphafold
example: 
    - parallel -j8 python app.py chain_af_summary ::: {0..240}
    - python app.py chain_af_summary 0
functions:
    - retrieve summary of alphafold prediction to table chain_af_summary
'''

from src.helper import *
from bioomics import ProcessJson

def parse_data(row, params):
    res = {}
    indir = os.path.join(params['predict_dir'], row['endpoint'])
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
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        SELECT a.relaxed_pdb
        FROM {table_name}     a
        LEFT JOIN pdb_chainid c ON a.chain_id = c.chain_id
        LEFT JOIN meta_pdb    p ON c.pdb_id = p.pdb_id
        WHERE p.pdb_group = {pdb_group}
            AND a.relaxed_pdb IS NOT NULL
    ;"""
    curr_rows = QueryComplex(params['verbose']).list_data(query)
    pdb_files = [i['relaxed_pdb'] for i in curr_rows]

    query = f"""
        SELECT a.pdb_id, a.chain_id, a.endpoint
        FROM chain_af      a
        LEFT JOIN meta_pdb p ON a.pdb_id = p.pdb_id
        WHERE p.pdb_group = {pdb_group}
            AND a.endpoint IS NOT NULL
    ;"""
    rows = QueryComplex(params['verbose']).list_data(query)
    meta['chains'] = len(rows)
    print(f"Filter {len(rows)} PDB data...")
    for row in rows:
        ranked = parse_data(row, params)
        for ranking, log in ranked.items():
            relaxed_pdb = log.get('relaxed_pdb')
            if relaxed_pdb and relaxed_pdb not in pdb_files:
                yield (
                    row['chain_id'],
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
        'chunk_size': 5,
        'table_name': 'chain_af_summary',
    })
    rows = retrieve_rows(params, meta)
    bc = BuildComplex(params['verbose'], params['chunk_size'])
    bc.insert_batch_records(rows, params['table_name'])

    footer(meta)

'''
requirements:
    - table chain_alphafold
example: 
    - parallel -j8 python app.py chain_afsum ::: {0..117}
    - python app.py chain_afsum 0
functions:
    - retrieve summary of alphafold prediction to table chain_afsum
'''

from src.helper import *
from bioomics import ProcessJson

def parse_data(row, args):
    res = {}
    indir = os.path.join(args['predict_dir'], row['endpoint'])
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
def retrieve_rows(args, meta):
    pdb_group = args['pdb_group']
    table_name = args['table_name']
    query = f"""
        SELECT a.relaxed_pdb
        FROM {table_name}     a
        LEFT JOIN pdb_chainid c ON a.chain_id = c.chain_id
        LEFT JOIN meta_pdb    p ON c.pdb_id = p.pdb_id
        WHERE p.pdb_group = {pdb_group}
            AND a.relative_pdb IS NOT NULL
    ;"""
    curr_rows = QueryComplex().list_data(query)
    pdb_files = [i['relaxed_pdb'] for i in curr_rows]

    query = f"""
        SELECT c.pdb_id, a.chain_id, a.endpoint
        FROM chain_alphafold a
        LEFT JOIN pdb_chainid c ON a.chain_id = c.chain_id
        LEFT JOIN meta_pdb p ON c.pdb_id = p.pdb_id
        WHERE p.pdb_group = {pdb_group}
            AND a.endpoint IS NOT NULL
    ;"""
    rows = QueryComplex().list_data(query)
    meta['chains'] = len(rows)
    print(f"Filter {len(rows)} PDB data...")
    for row in rows:
        ranked = parse_data(row, args)
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
    args = {
        'pdb_group': sys.argv[1],
        'verbose': True,
        'predict_dir': os.getenv('predict_dir'),
        'table_name': 'chain_afsum',
    }
    rows = retrieve_rows(args, meta)
    bc = BuildComplex(rows, args['verbose'])
    bc.insert_batch_records(args['table_name'])

    footer(meta)

'''
requirements:
    - chain_alphafold
example: 
    - parallel -j8 python app.py chain_tmalign ::: {0..117}
    - python app.py chain_tmalign 0
functions:
    - calculate TM-align score between model and sample structure
'''

from src.helper import *
from bioomics import EvalPred

def pull_raw_pdb(params, meta):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        SELECT c.pdb_id, c.chain_id, c.chain_pdb AS native_pdb
        FROM meta_chain_pdb c 
        LEFT JOIN meta_pdb  p ON c.pdb_id = p.pdb_id
        WHERE p.pdb_group = {pdb_group}
            AND c.chain_pdb IS NOT NULL
            AND c.chain_pdb NOT IN (
                SELECT native_pdb FROM {table_name}
            )
    ;"""
    rows = QueryComplex(params['verbose']).list_data(query)
    meta['rows_raw'] = len(rows)
    print(f"{len(rows)} PDB files ready...")
    return rows

def pull_fixed_pdb(params, meta):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        SELECT c.pdb_id, c.chain_id, c.fixed_pdb AS native_pdb
        FROM meta_chain_relaxed c 
        LEFT JOIN meta_pdb  p ON c.pdb_id = p.pdb_id
        WHERE p.pdb_group = {pdb_group}
            AND c.fixed_pdb IS NOT NULL
            AND c.fixed_pdb NOT IN (
                SELECT native_pdb FROM {table_name}
            )
    ;"""
    rows = QueryComplex(params['verbose']).list_data(query)
    meta['rows_fixed'] = len(rows)
    print(f"{len(rows)} PDB files ready...")
    return rows

def pull_relaxed_pdb(params, meta):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        SELECT c.pdb_id, c.chain_id, c.relaxed_pdb AS native_pdb
        FROM meta_chain_relaxed c 
        LEFT JOIN meta_pdb  p ON c.pdb_id = p.pdb_id
        WHERE p.pdb_group = {pdb_group}
            AND c.relaxed_pdb IS NOT NULL
            AND c.relaxed_pdb NOT IN (
                SELECT native_pdb FROM {table_name}
            )
    ;"""
    rows = QueryComplex(params['verbose']).list_data(query)
    meta['rows_relaxed'] = len(rows)
    print(f"{len(rows)} PDB files ready...")
    return rows

def calculate_tmalign(rows:list, pred:dict, params, meta):
    ep = EvalPred(params['mmalign_dir'])

    for row in rows:
        chain_id = row['chain_id']
        native_pdb = os.path.join(params['data_dir'], row['native_pdb'])
        pdb_status = 'Raw'
        if row['native_pdb'].endswith('fixed.pdb'):
            pdb_status = 'Fixed'
        elif row['native_pdb'].endswith('relaxed.pdb'):
            pdb_status = 'Relaxed'

        if os.path.isfile(native_pdb) and chain_id in pred:
            pred_pdbs = pred[chain_id]['relative_pdb']
            for pred_pdb in pred_pdbs.tolist():
                model_pdb = os.path.join(params['predict_dir'], pred_pdb)
                res = {}
                if os.path.isfile(native_pdb) and 'rank_001' in pred_pdb:
                    res = ep.TMalign(native_pdb, model_pdb)
                    if res:
                        meta['succeed'] += 1
                    else:
                        meta['fail'] += 1
                    rec = (row['pdb_id'], chain_id, res.get('TM1'), res.get('TM2'), \
                        res.get('RMSD'), pdb_status, row['native_pdb'], pred_pdb)
                    if params['verbose']:
                        print(rec)
                    yield rec


if __name__ == "__main__":
    params = {
        'pdb_group': sys.argv[1],
        'verbose': True,
        'chunk_size': 50,
        'data_dir': os.getenv('data_dir'),
        'predict_dir': os.getenv('predict_dir'),
        'mmalign_dir': os.getenv('mmalign_dir'),
        'table_name': 'chain_tmalign',
    }
    pred = QueryComplex(params['verbose']).parse_chain_prediction()
    pool = [pull_raw_pdb, pull_fixed_pdb, pull_relaxed_pdb,]
    for func in pool:
        rows = func(params, meta)
        record_iter = calculate_tmalign(rows, pred, params, meta)
        bc = BuildComplex(record_iter, params['verbose'], params['chunk_size'])
        bc.insert_batch_records(params['table_name'])

    footer(meta)
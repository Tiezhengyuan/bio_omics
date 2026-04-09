'''
requirements:
    - chain_plddt
example: 
    - parallel -j8 python app.py chain_rmsd ::: {0..237}
    - python app.py chain_rmsd 0
functions:
    - similarity between model and sample structure
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
    print(f"Detect {len(rows)} PDB files...")
    return rows

def pull_fixed_pdb(params, meta):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        SELECT c.pdb_id, c.chain_id, c.fixed_pdb AS native_pdb
        FROM meta_chain_relaxed c 
        LEFT JOIN meta_pdb      p ON c.pdb_id = p.pdb_id
        WHERE p.pdb_group = {pdb_group}
            AND c.fixed_pdb IS NOT NULL
            AND c.fixed_pdb NOT IN (
                SELECT native_pdb FROM {table_name}
            )
    ;"""
    rows = QueryComplex(params['verbose']).list_data(query)
    meta['rows_fixed'] = len(rows)
    print(f"Detect {len(rows)} new fixed PDB files...")
    return rows

def pull_relaxed_pdb(params, meta):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        SELECT c.pdb_id, c.chain_id, c.relaxed_pdb AS native_pdb
        FROM meta_chain_relaxed c 
        LEFT JOIN meta_pdb      p ON c.pdb_id = p.pdb_id
        WHERE p.pdb_group = {pdb_group}
            AND c.relaxed_pdb IS NOT NULL
            AND c.relaxed_pdb NOT IN (
                SELECT native_pdb FROM {table_name}
            )
    ;"""
    rows = QueryComplex(params['verbose']).list_data(query)
    meta['rows_relaxed'] = len(rows)
    print(f"Detect {len(rows)} new relaxed PDB files...")
    return rows

def calculate_rmsd(rows:list, pred:dict, params, meta):
    am = EvalPred()
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
                rmsd = None
                if os.path.isfile(native_pdb) and 'rank_001' in pred_pdb:
                    try:
                        rmsd = am.rmsd(native_pdb, model_pdb)
                        if rmsd:
                            meta['succeed'] += 1
                        else:
                            meta['null'] += 1
                    except Exception as e:
                        meta['fail'] += 1
                yield (row['pdb_id'], chain_id, rmsd, pdb_status, row['native_pdb'], pred_pdb)



if __name__ == "__main__":
    params.update({
        'chunk_size': 50,
        'table_name': 'chain_rmsd',
    })
    pool = [
        pull_raw_pdb(params, meta),
        pull_fixed_pdb(params, meta),
        pull_relaxed_pdb(params, meta),
    ]
    pred = QueryComplex(params['verbose']).parse_chain_prediction()
    for rows in pool:
        record_iter = calculate_rmsd(rows, pred, params, meta)
        bc = BuildComplex(record_iter, params['verbose'], params['chunk_size'])
        bc.insert_batch_records(params['table_name'])

    footer(meta)
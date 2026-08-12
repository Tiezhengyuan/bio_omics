'''
functions:
    - used by 8_summary.py
'''

from src.helper import *
from bioomics import EvalPred

def retrieve_pred(params, meta):
    '''
    pdb_id, chain_id, endpoint
    '''
    table_name = params['predict_table_name']
    col_name = params['input_col_name']
    pdb_group = params['pdb_group']
    query = f"""
        SELECT a.pdb_id, a.{col_name}, a.endpoint
        FROM {table_name}  a
        LEFT JOIN meta_pdb p ON a.pdb_id = p.pdb_id
        WHERE p.pdb_group = {pdb_group}
            AND a.endpoint IS NOT NULL
    ;"""
    _rows = QueryComplex(params['verbose']).list_data(query)
    rows, del_rows = [], []
    for row in _rows:
        outdir = params['predict_dir'] / row['endpoint']
        if outdir.is_dir():
            row['outdir'] = outdir
            rows.append(row)
        else:
            del_rows.append(row)
    meta[params['category']] = len(rows)
    print(f"Detect {len(rows)} predictions by {params['category']}...")
    return rows

def pred_summary(params):
    params_table = params['params_table_name']
    table_name = params['table_name']
    pdb_group = params['pdb_group']
    query = f"""
        SELECT a.relaxed_pdb_name
        FROM {table_name}        a
        LEFT JOIN {params_table} b ON a.input_id = b.input_id
        LEFT JOIN meta_pdb       p ON b.pdb_id = p.pdb_id
        WHERE p.pdb_group = {pdb_group}
            AND a.relaxed_pdb_name IS NOT NULL
    ;"""
    curr_rows = QueryComplex(params['verbose']).list_data(query)
    relaxed_pdb_names = [i['relaxed_pdb_name'] for i in curr_rows]
    return relaxed_pdb_names


def chain_raw_model_pdb(params):
    pdb_status = 'raw'
    pdb_group = params['pdb_group']
    summary = params['summary_table_name']
    table_name = params['table_name']
    query = f"""
        SELECT c.pdb_id, '{pdb_status}' AS pdb_status, 
            c.chain_id AS input_id, s.category,
            c.relative_pdb AS native_pdb, s.relaxed_pdb AS model_pdb
        FROM chain_pdb       c 
        LEFT JOIN meta_pdb   p ON c.pdb_id = p.pdb_id
        LEFT JOIN {summary} s ON c.chain_id = s.input_id
        WHERE p.pdb_group = {pdb_group}
            AND s.ranking = 1
            AND c.relative_pdb IS NOT NULL
            AND c.relative_pdb NOT IN (
                SELECT native_pdb FROM {table_name}
            )
    ;"""
    return pdb_status, 'chain', query

def chain_fixed_model_pdb(params):
    pdb_status = 'fixed'
    pdb_group = params['pdb_group']
    summary = params['summary_table_name']
    table_name = params['table_name']
    query = f"""
        SELECT c.pdb_id, '{pdb_status}' AS pdb_status, 
            c.chain_id AS input_id, s.category,
            c.fixed_pdb AS native_pdb, s.relaxed_pdb AS model_pdb
        FROM chain_fixed    c 
        LEFT JOIN meta_pdb  p ON c.pdb_id = p.pdb_id
        LEFT JOIN {summary} s ON c.chain_id = s.input_id
        WHERE p.pdb_group = {pdb_group}
            AND s.ranking = 1
            AND c.fixed_pdb IS NOT NULL
            AND c.fixed_pdb NOT IN (
                SELECT native_pdb FROM {table_name}
            )
    ;"""
    return pdb_status, 'chain', query

def chain_relaxed_model_pdb(params):
    pdb_status = 'relaxed'
    pdb_group = params['pdb_group']
    summary = params['summary_table_name']
    table_name = params['table_name']
    query = f"""
        SELECT c.pdb_id, '{pdb_status}' AS pdb_status, 
            c.chain_id AS input_id, s.category,
            c.relaxed_pdb AS native_pdb, s.relaxed_pdb AS model_pdb
        FROM chain_relaxed  c 
        LEFT JOIN meta_pdb  p ON c.pdb_id = p.pdb_id
        LEFT JOIN {summary} s ON c.chain_id = s.input_id
        WHERE p.pdb_group = {pdb_group}
            AND s.ranking = 1
            AND c.relaxed_pdb IS NOT NULL
            AND c.relaxed_pdb NOT IN (
                SELECT native_pdb FROM {table_name}
            )
    ;"""
    return pdb_status, 'chain', query

def combo_raw_model_pdb(params, combo):
    pdb_status = 'raw'
    pdb_group = params['pdb_group']
    summary = params['summary_table_name']
    table_name = params['table_name']
    query = f"""
        SELECT c.pdb_id, '{pdb_status}' AS pdb_status,
            c.combo_id AS input_id, s.category,
            c.combo_pdb AS native_pdb, s.relaxed_pdb AS model_pdb
        FROM {combo}       c 
        LEFT JOIN meta_pdb   p ON c.pdb_id = p.pdb_id
        LEFT JOIN {summary} s ON c.combo_id = s.input_id
        WHERE p.pdb_group = {pdb_group}
            AND s.ranking = 1
            AND c.combo_pdb IS NOT NULL
            AND c.combo_pdb NOT IN (
                SELECT native_pdb FROM {table_name}
            )
    ;"""
    return pdb_status, combo, query

def combo_fixed_model_pdb(params, combo):
    pdb_status = 'fixed'
    pdb_group = params['pdb_group']
    summary = params['summary_table_name']
    table_name = params['table_name']
    query = f"""
        SELECT c.pdb_id, '{pdb_status}' AS pdb_status,
            c.combo_id AS input_id, s.category,
            c.fixed_pdb AS native_pdb, s.relaxed_pdb AS model_pdb
        FROM {combo}_fixed  c 
        LEFT JOIN meta_pdb  p ON c.pdb_id = p.pdb_id
        LEFT JOIN {summary} s ON c.combo_id = s.input_id
        WHERE p.pdb_group = {pdb_group}
            AND s.ranking = 1
            AND c.fixed_pdb IS NOT NULL
            AND c.fixed_pdb NOT IN (
                SELECT native_pdb FROM {table_name}
            )
    ;"""
    return pdb_status, combo, query
 
def combo_relaxed_model_pdb(params, combo):
    pdb_status = 'relaxed'
    pdb_group = params['pdb_group']
    summary = params['summary_table_name']
    table_name = params['table_name']
    query = f"""
        SELECT c.pdb_id, '{pdb_status}' AS pdb_status,
            c.combo_id AS input_id, s.category,
            c.relaxed_pdb AS native_pdb, s.relaxed_pdb AS model_pdb
        FROM {combo}_relaxed c 
        LEFT JOIN meta_pdb   p ON c.pdb_id = p.pdb_id
        LEFT JOIN {summary}  s ON c.combo_id = s.input_id
        WHERE p.pdb_group = {pdb_group}
            AND s.ranking = 1
            AND c.relaxed_pdb IS NOT NULL
            AND c.relaxed_pdb NOT IN (
                SELECT native_pdb FROM {table_name}
            )
    ;"""
    return pdb_status, combo, query

def calculate_rmsd(rows:list, params, meta):
    prefix = f"{params['pdb_status']}_{params['category']}"
    am = EvalPred()
    
    for row in rows:
        native_pdb = params['simulate_dir'] / row['native_pdb'] \
            if params['pdb_status']=='relaxed' \
            else params['output_pdb_dir'] / row['native_pdb']
        model_pdb = params['predict_dir'] / row['model_pdb']
        if native_pdb.is_file() and model_pdb.is_file():
            try:
                rmsd = am.rmsd(native_pdb, model_pdb)
                if rmsd:
                    meta[f'{prefix}_succeed'] += 1
                else:
                    meta[f'{prefix}_null'] += 1
                # build record
                rec = (
                    row['pdb_id'],
                    row['category'],
                    row['input_id'],
                    row['pdb_status'],
                    rmsd,
                    row['native_pdb'],
                    row['model_pdb'],
                )
                yield rec
            except Exception as e:
                meta[f'{prefix}_fail'] += 1
                pass

def calculate_tmalign(rows:list, params, meta):
    prefix = f"{params['pdb_status']}_{params['category']}"
    ep = EvalPred(params['mmalign_dir'])
    
    for row in rows:
        native_pdb = params['simulate_dir'] / row['native_pdb'] \
            if params['pdb_status']=='relaxed' \
            else params['output_pdb_dir'] / row['native_pdb']
        model_pdb = params['predict_dir'] / row['model_pdb']
        if native_pdb.is_file() and model_pdb.is_file():
            try:
                res = ep.TMalign(native_pdb, model_pdb)
                if res:
                    meta[f'{prefix}_succeed'] += 1
                else:
                    meta[f'{[prefix]}_fail'] += 1
                rec = (
                    row['pdb_id'],
                    row['category'],
                    row['input_id'],
                    row['pdb_status'],
                    res.get('RMSD'),
                    res.get('TM1'),
                    res.get('TM2'),
                    row['native_pdb'],
                    row['model_pdb'],
                )
                yield rec
            except Exception as e:
                meta[f'{prefix}_fail'] += 1
                pass
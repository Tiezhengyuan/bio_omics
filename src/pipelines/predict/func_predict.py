'''
fuctions for structure prediction
'''
import asyncio
from src.helper import *
from bioomics import Predictor, WriteTmp

# def validate_predict_data(params, meta):
#     col_name = 'endpoint'
#     pdb_group = params['pdb_group']
#     table_name = params['table_name']
#     # scan no existing files
#     print(f"Scan {table_name}.{col_name}, group={pdb_group}...")
#     query = f"""
#         SELECT S.pdb_id, S.{col_name}
#         FROM {table_name}  S
#         LEFT JOIN meta_pdb P ON S.pdb_id = P.pdb_id
#         WHERE P.pdb_group = {pdb_group}
#             AND S.{col_name} IS NOT NULL
#     ;"""
#     rows = QueryComplex(params['verbose']).list_data(query)
#     del_pool = []
#     for row in rows:
#         outdir = params['predict_dir'] / row[col_name]
#         if not outdir.is_dir():
#             print(f"NOT EXISTING: {outdir}.")
#             del_pool.append(row[col_name])
#     # delete rows
#     dc = DeleteComplex(params['verbose'])
#     num_del = dc.remove_null_rows(table_name, col_name, del_pool)
#     meta['delete_rows'] = num_del

def validate_file_rows(col_name, params, meta):
    params_table = params['params_table_name']
    table_name = params['table_name']
    pdb_group = params['pdb_group']
    query = f"""
        SELECT a.input_id, b.outdir, a.{col_name}
        FROM {table_name}        a
        LEFT JOIN {params_table} b ON a.input_id = b.input_id
        LEFT JOIN meta_pdb       m ON b.pdb_id = m.pdb_id
        WHERE m.pdb_group = {pdb_group}
    """
    rows = QueryComplex(params['verbose']).list_data(query)
    del_rows = []
    for row in rows:
        path = Path(row['outdir']) / row[col_name]
        if not path.is_file():
            del_rows.append(row['input_id'])
            
    # delete null rows
    if del_rows:
        print(f"Delete {len(del_rows)} *.asm...")
        dc = DeleteComplex(params['verbose'])
        dc.remove_null_rows(table_name, 'input_id', del_rows)
        meta['delete_rows'] = len(del_rows)

        
def retrieve_chain(params, outdir_name:str):
    input_table_name = 'label_chain'
    input_col_name = 'chain_id'
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        WITH tmp AS(
            SELECT DISTINCT a.first_chain_id AS input_id
            FROM {input_table_name} a
            LEFT JOIN meta_pdb      m ON a.pdb_id = m.pdb_id
            WHERE m.pdb_group = {pdb_group}
                AND NOT EXISTS (
                    SELECT 1 FROM {table_name} c
                    WHERE c.input_id = a.first_chain_id
                )
        )
        SELECT t.input_id, f.pdb_id, f.relative_faa
        FROM tmp t, pdb_chain_seq f
        WHERE t.input_id = f.chain_id
    """
    cat_params = dict(params)
    cat_params.update({
        'input_table_name': input_table_name,
        'input_col_name': input_col_name,
        'input_dir': params['output_pdb_dir'],
        'category': 'chain',
        'outdir_name': outdir_name,
        'outdir': params['predict_dir'] / outdir_name,
    })
    return cat_params, query

def retrieve_combo(params, outdir_name:str, combo:str):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    input_table_name = 'label_complex'
    input_col_name = 'combo_id'
    query = f"""
        SELECT DISTINCT a.combo_id AS input_id, a.pdb_id,
            c.combo_fasta AS relative_faa
        FROM label_complex a
        LEFT JOIN {combo}  c ON a.combo_id = c.combo_id
        LEFT JOIN meta_pdb p ON a.pdb_id = p.pdb_id
        WHERE p.pdb_group = {pdb_group}
            AND c.combo_fasta IS NOT NULL
            AND NOT EXISTS (
                SELECT 1 FROM {table_name} b
                WHERE a.combo_id = b.input_id
            )
    """
    cat_params = dict(params)
    cat_params.update({
        'input_table_name': input_table_name,
        'input_col_name': input_col_name,
        'input_dir': params['output_pdb_dir'],
        'outdir_name': outdir_name,
        'category': combo,
        'outdir_name': outdir_name,
        'outdir': params['predict_dir'] /outdir_name,
    })
    return cat_params, query

def build_curr_params(row, params):
    # update output_predict_dir for exporting outputs
    pdb_id, input_id = row['pdb_id'], row['input_id']
    endpoint = Path(params['outdir_name']) / pdb_id[:2] / pdb_id / input_id
    
    curr_params = dict(params)
    curr_params.update({
        'predicted_at': datetime.now(),
        'host': params['host'],
        'pdb_id': pdb_id,
        'input_id': input_id,
        'endpoint': endpoint,
        'faa_file': params['input_dir'] / row['relative_faa'],
        'outdir': params['predict_dir'] / endpoint,
    })
    if params['verbose']:
        print('Arguments for prediction: ', curr_params)
    return curr_params

async def run_colabfold_search(input):
    row, params, meta = input
    # check if MSA is available
    if params['overwrite'] or not row['a3m_file'].is_file():
        outdir = row['outdir']
        tmp_name = f"colabfold_search_{row['input_id']}"
        with WriteTmp(params['tmp_dir'], outdir, tmp_name) as writer:
            # run MSA    
            search_params = {
                'tools_dir': params['tools_dir'],
                'query_faa': row['query_faa'],
                'outdir': writer.temp_dir_path,
                'colabfold_db': params['colabfold_db'],
            }
            res = await asyncio.to_thread(
                Predictor.colabfold_search,
                search_params,
            )
            if res:
                meta['msa_good'] += 1
            else:
                meta['msa_fail'] += 1
    else:
        meta['no_a3m'] += 1
    return row, params, meta
    

async def run_colabfold_msa(input):
    row, params, meta = input
    # check if MSA is available
    if params['overwrite'] or not row['a3m_file'].is_file():
        outdir = row['outdir']
        tmp_name = f"colabfold_msa_{row['input_id']}"
        with WriteTmp(params['tmp_dir'], outdir, tmp_name) as writer:
            # run MSA 
            api_params = {
                'tools_dir': params['tools_dir'],
                'query_faa': row['query_faa'],
                'outdir': writer.temp_dir_path,
            }
            res = await asyncio.to_thread(
                Predictor.colabfold_batch_msa,
                api_params
            )
            if res:
                meta['msa_good'] += 1
            else:
                meta['msa_fail'] += 1
    
        # delete some directories for saving space
        for dir_suffix in ('env', 'pairgreedy', '_unpaired'):
            _dir = Path(row['outdir'])  / f"{row['input_id']}_{dir_suffix}"
            if _dir.is_dir():
                shutil.rmtree(_dir)
    else:
        meta['no_a3m'] += 1
    return row, params, meta

def retrieve_msa(params, meta):
    params_table = params['params_table_name']
    msa_table = params['msa_table_name']
    pdb_group = params['pdb_group']
    predictor = params['predictor']
    query = f"""
        SELECT b.pdb_id, b.input_id, b.outdir, a.a3m_name
        FROM {params_table}   b 
        LEFT JOIn {msa_table} a ON b.input_id = a.input_id
        LEFT JOIN meta_pdb    p ON b.pdb_id = p.pdb_id
        WHERE p.pdb_group = {pdb_group}
            AND a.a3m_name IS NOT NULL
            AND b.predictor = '{predictor}'
        ORDER BY b.category
    ;"""
    rows = QueryComplex(params['verbose']).list_data(query)
    meta['msa_rows'] = len(rows)
    
    frows = []
    for row in rows:
        row['outdir'] = Path(row['outdir'])
        row['input_file'] = row['outdir'] / row['a3m_name']
        if row['input_file'].is_file():
            frows.append(row)
    print(f"Detect {len(frows)} *.asm...")

    return frows


def run_colabfold_batch(row, params, meta):
    # check if prediction is done
    done_name = f"{row['input_id']}.done.txt"
    done_file = row['outdir'] / done_name
    if not done_file.is_file():
        outdir = row['outdir']
        with WriteTmp(params['tmp_dir'], outdir, 'colabfold_batch') as writer:
            #keys: tools_dir, input_file, outdir
            tmp_outdir = writer.temp_dir_path
            p = Predictor(
                row['input_file'],
                tmp_outdir,
                params['tools_dir'],
                params['verbose']
            )
            res =  p.colabfold_batch(params)
            if res:
                meta['succeed'] += 1
            else:
                meta['fail'] += 1

            # delete *_env for saving space
            _dir = tmp_outdir  / f"{row['input_id']}_all"
            if _dir.is_dir():
                shutil.rmtree(_dir)
            

'''
fuctions for structure prediction
'''

from src.helper import *
from bioomics import Predictor, TABLE_COLUMNS, ProcessFasta, WriteTmp

def validate_predict_data(params, meta):
    col_name = 'endpoint'
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    # scan no existing files
    print(f"Scan {table_name}.{col_name}, group={pdb_group}...")
    query = f"""
        SELECT S.pdb_id, S.{col_name}
        FROM {table_name}  S
        LEFT JOIN meta_pdb P ON S.pdb_id = P.pdb_id
        WHERE P.pdb_group = {pdb_group}
            AND S.{col_name} IS NOT NULL
    ;"""
    rows = QueryComplex(params['verbose']).list_data(query)
    del_pool = []
    for row in rows:
        outdir = params['predict_dir'] / row[col_name]
        if not outdir.is_dir():
            print(f"NOT EXISTING: {outdir}.")
            del_pool.append(row[col_name])
    # delete rows
    dc = DeleteComplex(params['verbose'])
    num_del = dc.remove_null_rows(table_name, col_name, del_pool)
    meta['delete_rows'] = num_del
        
@header
def retrieve_chain(params, meta):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        SELECT DISTINCT F.pdb_id, F.chain_id, F.relative_faa
        FROM label_chain        S
        LEFT JOIN pdb_chain_seq F ON S.first_chain_id = F.chain_id
        LEFT JOIN meta_pdb      P ON S.pdb_id = P.pdb_id
        WHERE P.pdb_group = {pdb_group}
            AND F.relative_faa IS NOT NULL
            AND S.first_chain_id IS NOT NULL
            AND S.first_chain_id NOT IN (
                SELECT chain_id FROM {table_name}
                WHERE chain_id IS NOT NULL
            )
    """
    rows = QueryComplex(params['verbose']).list_data(query)
    print(f"{len(rows)} PDB files ready...")
    meta['rows'] += len(rows)
    return rows

def chain_params(row, params):
    pdb_id, chain_id = row['pdb_id'], row['chain_id']
    endpoint = Path(params['table_name']) / pdb_id[:2] / pdb_id / chain_id
    # update output_predict_dir for exporting outputs
    outdir = params['predict_dir'] / endpoint
    outdir.mkdir(parents=True, exist_ok=True)
    params.update({
        'predicted_at': datetime.now(),
        'host': params['host'],
        'pdb_id': pdb_id,
        'chain_id': chain_id,
        'endpoint': str(endpoint),
        'faa_file': params['output_pdb_dir'] / row['relative_faa'],
        'outdir': outdir,
    })
    params['input_id'] = params['chain_id']
    if params['verbose']:
        print('Arguments for prediction: ', params)

@header
def retrieve_combo(params, meta):
    qc = QueryComplex(params['verbose'])
    rows = qc.retrieve_combo_fasta(params['pdb_group'], params['combo'], params['table_name'])
    meta['rows'] += len(rows)

    # antibody chain should be included in combo
    ab_rows = []
    for row in rows:
        row['faa_file'] = os.path.join(params['output_pdb_dir'], row['combo_fasta'])
        # predict short sequence for saving time
        records = ProcessFasta(row['faa_file']).load_fasta()
        if isinstance(records, list):
            total_aa = np.sum([len(rec.seq) for rec in records])
            if total_aa < 1000:
                print(row['combo_id'], total_aa, end='\t')
                ab_rows.append(row)
    meta['predict_rows'] += len(ab_rows)
    print(f"Get {len(ab_rows)} PDB files for prediction.")
    return ab_rows

def combo_params(row, params):
    pdb_id, combo_id = row['pdb_id'], row['combo_id']
    endpoint = Path(params['table_name']) / pdb_id[:2] / pdb_id / combo_id
    outdir = params['predict_dir'] / endpoint
    outdir.mkdir(parents=True, exist_ok=True)
    params.update({
        'predicted_at': datetime.now(),
        'host': params['host'],
        'pdb_id': pdb_id,
        'combo_id': combo_id,
        'endpoint': str(endpoint), # string type for insertion
        'faa_file': params['output_pdb_dir'] / row['combo_fasta'],
        'outdir':outdir,
    })
    params['input_id'] = params['combo_id']
    if params['verbose']:
        print('Arguments for prediction: ', params)

def run_colabfold_batch(params, meta):
    # store data in RAM
    outdir = params['outdir']
    with WriteTmp(params['tmp_dir'], outdir, params['table_name']) as writer:
        # run predictor
        params['outdir'] = writer.temp_dir_path
        res = Predictor(params).colabfold_batch()

        if res:
            # check if prediction is done
            done_file = params['outdir'] / f"{params['input_id']}.done.txt"
            if done_file.is_file():
                meta['succeed'] += 1
            else:
                params['endpoint'] = None
                meta['fail'] += 1

            # insertion
            rec = [params.get(col) for col in TABLE_COLUMNS[params['table_name']]]
            if params['verbose']:
                print('### Try to insert ', rec)
            bc = BuildComplex(params['verbose'])
            bc.insert_batch_records([tuple(rec),], params['table_name'])

            # delete *_env for saving space
            if 'endpoint' in params:
                _env_dir = params['outdir']  / f"{params['input_id']}_env"
                if _env_dir.is_dir():
                    shutil.rmtree(_env_dir)
        else:
            print(f"WARNING: {params['tool_name']} failed!")
#!/usr/bin/python
'''
requirements:
    - table view_antibody, label_chain, lablel_complex
example:
    - python app.py predict_npr1
function:
    - design complex
'''
from src.helper import *
from bioomics import Predictor

def pool_data(params):
    target_name = params['target_name']
    query = f"""
        SELECT category, target_name, design_id, chain_labels, af_faa
        FROM design_dimer
        WHERE target_name = '{target_name}'
        UNION
        SELECT category, target_name, design_id, chain_labels, af_faa
        FROM design_multimer
        WHERE target_name = '{target_name}'
    ;"""
    rows = QueryComplex(params['verbose']).list_data(query)
    df = pd.DataFrame(rows).set_index('design_id', drop=False)
    return df

# status == stop
def retrieve_stopped(params, df0) -> list:
    table_name = params['table_name']
    query = f"""
        SELECT design_id FROM {table_name}
        WHERE host = '{params['host']}' AND status = 'stop'
    ;"""
    rows = QueryComplex(params['verbose']).list_data(query)
    ids = set([row['design_id'] for row in rows])
    df1 = df0.loc[list(ids.intersection(df0.index))]
    meta['unfinished _rows'] = len(df1)

    # update status
    res = []
    uc = UpdateComplex(params['table_name'], params['verbose'])
    for _, row in df1.iterrows():
        uc.design_status(params['predictor'], row['design_id'], 'run')
        res.append(row)
    return res

# 1. there are results but no record in table predict
# 2. new data
def detect_unrecorded(params, df0) -> list:
    # collect recorded data
    q = QueryComplex(params['verbose'])
    ids = q.unique_values(params['table_name'], 'design_id')

    # unrecorded
    unrecorded_ids = []
    output_dir = os.path.join(params['design_dir'], params['table_name'])
    if os.path.isdir(output_dir):
        for endpoint_name in os.listdir(output_dir):
            path = os.path.join(output_dir, endpoint_name)
            for design_id in os.listdir(path):
                if design_id not in ids:
                    unrecorded_ids.append(design_id)
        df1 = df0[df0['design_id'].isin(unrecorded_ids)]
        if not df1.empty:
            return df1
    
    # retrieve new data
    df1 = df0[~df0['design_id'].isin(ids)]
    return df1

def retrieve_new(params, df1) -> list:
    res = []
    if not df1.empty:
        g = df1.groupby(['category', 'target_name'])
        for _, sub in g:
            sub1 = sub.sample(n=params['batch_size'], replace=False)
            for _, row in sub1.iterrows():
                # insert record first
                rec = (params['predictor'], row['design_id'], params['host'], 'run')
                bc = BuildComplex([rec,], params['verbose'])
                bc.insert_batch_records(params['table_name'])
                res.append(row)
    return res

if __name__ == "__main__":
    params.update({
        'design_dir': os.getenv('design_dir'),
        'table_name': 'predict',
        'predictor': 'AF2',
        'batch_size': 2,
        'num_recycle': '3',
        'relax': True,
        'target_name': 'NPR1',
    })
    df0= pool_data(params)
    print("Firstly, try unfinised previously...")
    rows = retrieve_stopped(params, df0)
    if not rows:
        print(f"Try to detect unrecorded or newly data...")
        df1 = detect_unrecorded(params, df0)
        rows = retrieve_new(params, df1)

    uc = UpdateComplex(params['table_name'], params['verbose'])
    for row in rows:
        t0 = datetime.now()
        endpoint = os.path.join(params['table_name'], row['chain_labels'], row['design_id'])
        outdir = os.path.join(params['design_dir'], endpoint)
        Path(outdir).mkdir(parents=True, exist_ok=True)

        # predict
        _params = {
            'faa_file': row['af_faa'],
            'output_dir': outdir,
        }
        _params.update(params)
        if params['verbose']:
            print('Arguments for prediction: ', _params)
        Predictor(_params).colabfold()

        # update record
        done_file = os.path.join(outdir, f"{row['design_id']}.done.txt")
        if os.path.isfile(done_file):
            uc.design_end_time(params['predictor'], row['design_id'], endpoint)
            meta['succeed'] += 1
            # delete *_env for saving space
            _env_dir = os.path.join(outdir, f"{row['design_id']}_env")
            if os.path.isdir(_env_dir):
                shutil.rmtree(_env_dir)
            #copy results to local server if that is remote
            copy_data(outdir)
        else:
            uc.design_status(params['predictor'], row['design_id'], 'failed')
            meta['fail'] += 1


    footer(params, meta)
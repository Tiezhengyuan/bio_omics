#!/usr/bin/python
'''
requirements:
    - table label_chain
example:
    - parallel -j8 python app.py chain_af ::: {0..240}
    - python app.py chain_af 0
function:
    - prefix monomer
    - insert data into table chain_af
'''

from src.helper import *
from bioomics import Predictor, TABLE_COLUMNS

@header
def retrieve_pdb(params, meta):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        SELECT DISTINCT F.pdb_id, F.chain_id, F.relative_path
        FROM label_chain    S
        LEFT JOIN chain_faa F ON S.first_chain_id = F.chain_id
        LEFT JOIN meta_pdb  P ON S.pdb_id = P.pdb_id
        WHERE P.pdb_group = {pdb_group}
            AND F.relative_path IS NOT NULL
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

if __name__ == "__main__":
    params.update({
        'table_name': 'chain_af',
        'num_recycle': '3',
        'relax': True,
        'max_iterations': '2000',
    })
    add_params(params)

    rows = retrieve_pdb(params, meta)
    for row in rows:
        pdb_id, chain_id = row['pdb_id'], row['chain_id']
        item = {
            'predicted_at': datetime.now(),
            'host': params['host'],
            'pdb_id': pdb_id,
            'chain_id': chain_id,
        }
        endpoint = os.path.join(params['table_name'], pdb_id[:2], pdb_id, chain_id)
        outdir = os.path.join(params['predict_dir'], endpoint)
        Path(outdir).mkdir(parents=True, exist_ok=True)

        # predict
        _params = {
            'faa_file': os.path.join(params['output_pdb_dir'], row['relative_path']),
            'output_dir': outdir,
        }
        _params.update(params)
        if params['verbose']:
            print('Arguments for prediction: ', _params)
        Predictor(_params).colabfold()

        # check if prediction is done
        done_file = os.path.join(outdir, f"{chain_id}.done.txt")
        if os.path.isfile(done_file):
            item['endpoint'] = endpoint
            meta['succeed'] += 1
        else:
            meta['fail'] += 1

        # insertion
        rec = [item.get(col) for col in TABLE_COLUMNS[params['table_name']]]
        if params['verbose']:
            print('### Try to insert ', rec)
        bc = BuildComplex(params['verbose'])
        bc.insert_batch_records([tuple(rec),], params['table_name'])

        # delete *_env for saving space
        if 'endpoint' in item:
            _env_dir = os.path.join(outdir, f"{chain_id}_env")
            if os.path.isdir(_env_dir):
                shutil.rmtree(_env_dir)
            #copy results to local server if that is remote
            copy_data(outdir)

    footer(meta, params)
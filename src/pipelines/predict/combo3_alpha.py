#!/usr/bin/python
'''
requirements:
    - meta_combo3
example:
    - parallel --delay 5 -j 8 python app.py combo3_alpha ::: {0..99}
    - python app.py combo3_alpha
function:
    - predict dimer complex
    - insert data into table combo3_alpha
'''


from src.helper import *
from bioomics import Predictor

@header
def retrieve_pdb(params, meta):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        SELECT b.pdb_id, b.combo_id, c.combo_fasta
        FROM pdb_complex       b 
        LEFT JOIN meta_combo3  c ON b.combo_id = c.combo_id 
        LEFT JOIN meta_pdb     p ON b.pdb_id = p.pdb_id
        WHERE b.category = 'HLAg'
            AND p.pdb_group = {pdb_group}
            AND b.combo_id NOT IN (
                SELECT combo_id FROM {table_name}
            )
    ;"""
    qc = QueryComplex(params['verbose'])
    rows = qc.list_data(query)
    print(f"{len(rows)} PDB files ready...")
    meta['rows'] += len(rows)

    # antibody chain should be included in combo
    ab_rows = []
    for row in rows:
        row['faa_file'] = os.path.join(params['data_dir'], row['combo_fasta'])
        # predict short sequence for saving time
        records = ProcessFasta(row['faa_file']).load_fasta()
        total_aa = np.sum([len(rec.seq) for rec in records])
        if total_aa < 1000:
            print(row['combo_id'], total_aa, end='\t')
            ab_rows.append(row)
    meta['predict_rows'] += len(ab_rows)
    print(f"Get {len(ab_rows)} PDB files for prediction.")
    return ab_rows


if __name__ == "__main__":
    params = {
        'verbose': True,
        'data_dir': os.getenv('data_dir'),
        'predict_dir': os.getenv('predict_dir'),
        'table_name': 'combo3_alphafold',
        'num_recycle': '3',
        'relax': True,
    }
    add_params(params)

    rows = retrieve_pdb(params, meta)
    host = socket.gethostname()
    for row in rows:
        t0 = datetime.now()
        pdb_id, combo_id = row['pdb_id'], row['combo_id']
        endpoint = os.path.join(params['table_name'], pdb_id[:2], pdb_id, combo_id)
        outdir = os.path.join(params['predict_dir'], endpoint)
        Dir(outdir).init_dir()

        # predict
        _params = {
            'faa_file': row['faa_file'],
            'output_dir': outdir,
        }
        _params.update(params)
        if params['verbose']:
            print('Arguments for prediction: ', _params)
        Predictor(_params).colabfold()

        # build record
        rec = [pdb_id, combo_id, None, host, t0]
        done_file = os.path.join(outdir, f"{combo_id}.done.txt")
        if os.path.isfile(done_file):
            meta['succeed'] += 1
            rec[2] = endpoint
        else:
            meta['fail'] += 1
        
        # insertion
        if params['verbose']:
            print('### Try to insert ', rec)
        bc = BuildComplex([tuple(rec),], params['verbose'])
        bc.insert_batch_records(params['table_name'])

        if rec[2]:
            # delete *_env for saving space
            _env_dir = os.path.join(outdir, f"{combo_id}_env")
            if os.path.isdir(_env_dir):
                shutil.rmtree(_env_dir)
            #copy results to local server if that is remote
            copy_data(outdir)
    footer(meta, params)
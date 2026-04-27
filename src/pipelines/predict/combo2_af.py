#!/usr/bin/python
'''
requirements:
    - combo2
example:
    - parallel -j10 python app.py combo2_af ::: {0..237}
    - python app.py combo2_af 0
    - python app.py combo2_af
function:
    - predict dimer complex
    - insert data into table combo2_alphafold
'''

from src.helper import *
from bioomics import Predictor, ProcessFasta

@header
def retrieve_pdb(params, meta):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        SELECT DISTINCT c.pdb_id, c.combo_id, c.combo_fasta
        FROM combo2        c
        LEFT JOIN meta_pdb p ON c.pdb_id = p.pdb_id
        WHERE p.pdb_group = {pdb_group}
            AND NOT EXISTS (
                SELECT 1 FROM {table_name} b
                WHERE b.combo_id = c.combo_id
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
        if isinstance(records, list):
            total_aa = np.sum([len(rec.seq) for rec in records])
            if total_aa < 1000:
                print(row['combo_id'], total_aa, end='\t')
                ab_rows.append(row)
    meta['predict_rows'] += len(ab_rows)
    print(f"Get {len(ab_rows)} PDB files for prediction.")
    return ab_rows


if __name__ == "__main__":
    params.update({
        'table_name': 'combo2_alphafold',
        'num_recycle': '3',
        'relax': True,
    })
    add_params(params)

    rows = retrieve_pdb(params, meta)
    for row in rows:
        t0 = datetime.now()
        pdb_id, combo_id = row['pdb_id'], row['combo_id']
        endpoint = os.path.join(params['table_name'], pdb_id[:2], pdb_id, combo_id)
        outdir = os.path.join(params['predict_dir'], endpoint)
        Path(outdir).mkdir(parents=True, exist_ok=True)

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
        
        done_file = os.path.join(outdir, f"{combo_id}.done.txt")
        if os.path.isfile(done_file):
            meta['succeed'] += 1
        else:
            meta['fail'] += 1
            endpoint = None
        rec = [pdb_id, combo_id, endpoint, params['host'], t0]
        # insertion
        if params['verbose']:
            print('### Try to insert ', rec)
        bc = BuildComplex([tuple(rec),], params['verbose'])
        bc.insert_batch_records(params['table_name'])

        if endpoint:
            # delete *_env for saving space
            _env_dir = os.path.join(outdir, f"{combo_id}_env")
            if os.path.isdir(_env_dir):
                shutil.rmtree(_env_dir)
            #copy results to local server if that is remote
            copy_data(outdir)

    footer(meta, params)
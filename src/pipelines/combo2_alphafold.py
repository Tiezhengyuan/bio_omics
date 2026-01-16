#!/usr/bin/python
'''
requirements:
    - meta_combo2
example:
    - parallel -j8 python src/pipelines/combo2_alphafold.py ::: {0..117}
    - python src/pipelines/combo2_alphafold.py 0
function:
    - predict dimer complex
    - insert data into table combo2_alphafold
'''

from helper import *
from bioomics import Predictor

@header
def retrieve_pdb(args, meta):
    pdb_group = args['pdb_group']
    table_name = args['table_name']
    query = f"""
        SELECT C.pdb_id, C.combo_id, C.combo_fasta
        FROM meta_combo2   C
        LEFT JOIN meta_pdb P ON C.pdb_id=P.pdb_id
        WHERE P.pdb_group = {pdb_group}
            AND C.pdb_id IN (
                SELECT pdb_id from view_antibody
            )
            AND C.combo_id NOT IN (
                SELECT combo_id FROM {table_name}
            )
    ;"""
    qc = QueryComplex(args['verbose'])
    rows = qc.list_data(query)
    print(f"{len(rows)} PDB files ready...")
    meta['rows'] += len(rows)

    # antibody chain should be included in combo
    ab_rows = []
    ab_chains = qc.ab_chain_no()
    for row in rows:
        pdb_id = row['pdb_id']
        if pdb_id in ab_chains:
            combo_ids = set(row['combo_id'].split('-'))
            if combo_ids.intersection(ab_chains[pdb_id]):
                ab_rows.append(row)
    meta['antibody_rows'] += len(ab_rows)
    print(f"Filter {len(ab_rows)} PDB files including antibody chain.")
    return ab_rows


if __name__ == "__main__":
    args = {
        'pdb_group': sys.argv[1],
        'verbose': True,
        'data_dir': os.getenv('data_dir'),
        'predict_dir': os.getenv('predict_dir'),
        'table_name': 'combo2_alphafold',
        'num_recycle': '3',
        'relax': True,
        'max_iterations': '5000',
    }

    rows = retrieve_pdb(args, meta)
    for row in rows:
        t0 = datetime.now()
        pdb_id, combo_id = row['pdb_id'], row['combo_id']
        endpoint = os.path.join(args['table_name'], pdb_id[:2], pdb_id, combo_id)
        outdir = os.path.join(args['predict_dir'], endpoint)
        Dir(outdir).init_dir()

        # predict
        _args = {
            'faa_file': os.path.join(args['data_dir'], row['combo_fasta']),
            'output_dir': outdir,
        }
        _args.update(args)
        if args['verbose']:
            print('Arguments for prediction: ', _args)
        Predictor(_args).colabfold()

        # build record
        rec = [combo_id, None, t0]
        done_file = os.path.join(outdir, f"{combo_id}.done.txt")
        if os.path.isfile(done_file):
            meta['succeed'] += 1
            rec[1] = endpoint
        else:
            meta['fail'] += 1
        # insertion
        if args['verbose']:
            print('### Try to insert ', rec)
        bc = BuildComplex([tuple(rec),], args['verbose'])
        bc.insert_batch_records(args['table_name'])

    footer(meta)
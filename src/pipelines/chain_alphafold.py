#!/usr/bin/python
'''
requirements:
    - pdb fixed: table fixed_chain
example:
    - parallel -j4 python src/pipelines/chain_alphafold.py ::: {0..116}
    - python src/pipelines/chain_alphafold.py 0
function:
    - relax structure based on fixed pdb files
    - insert data into table chain_alphafold
'''
from collections import defaultdict
from datetime import datetime
from dotenv import load_dotenv
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import Dir, BuildComplex, QueryComplex, Predictor


def retrieve_pdb(args, meta):
    pdb_group = args['pdb_group']
    table_name = args['table_name']
    query = f"""
        SELECT P.pdb_id, F.chain_id, F.faa_file AS endpoint
        FROM meta_chain_faa F
        LEFT JOIN pdb_chainid I
        ON F.chain_id = I.chain_id
        LEFT JOIN meta_pdb P
        ON I.pdb_id = P.pdb_id
        LEFT JOIN pdb_chain_seq S
        ON I.chain_id = S.chain_id
        WHERE P.pdb_group = {pdb_group}
            AND F.faa_file IS NOT NULL
            AND S.chain_seq IS NOT NULL
            AND I.model_no = 0
            AND F.chain_id NOT IN (
                SELECT chain_id FROM {table_name}
            )
            AND P.pdb_id IN (
                SELECT pdb_id from view_antibody
            )
    ;"""
    rows = QueryComplex(args['verbose']).list_data(query)
    print(f"{len(rows)} PDB files ready...")
    meta['pdb_group'] = int(pdb_group)
    meta['rows'] += len(rows)
    return rows

if __name__ == "__main__":
    meta = defaultdict(int)
    load_dotenv()
    args = {
        'verbose': True,
        'data_dir': os.getenv('data_dir'),
        'predict_dir': os.getenv('predict_dir'),
        'table_name': 'chain_alphafold',
        'max_iterations': '5000',
        'num_recycle': '3',
    }
    args['pdb_group'] = sys.argv[1]
    rows = retrieve_pdb(args, meta)

    for row in rows:
        t0 = datetime.now()
        pdb_id, chain_id, = row['pdb_id'], row['chain_id'], 
        endpoint = os.path.join(args['table_name'], pdb_id[:2], pdb_id, chain_id)
        outdir = os.path.join(args['predict_dir'], endpoint)
        Dir(outdir).init_dir()

        # predict
        _args = {
            'faa_file': os.path.join(args['data_dir'], row['endpoint']),
            'output_dir': outdir,
        }
        _args.update(args)
        Predictor(_args).colabfold()

        # build record
        rec = [chain_id, None, t0]
        done_file = os.path.join(outdir, f"{chain_id}.done.txt")
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
    print(meta)
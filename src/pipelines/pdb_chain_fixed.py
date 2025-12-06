#!/usr/bin/python
'''
requirements:
    - table pdb_chainid
example:
    - parallel -j1 python src/pipelines/pdb_chain_fixed.py ::: {0..116}
    - python src/pipelines/pdb_chain_fixed.py 0
function:
    - fix pdb files
    - insert data into table pdb_chain_fixed
'''
from collections import defaultdict
from dotenv import load_dotenv
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics.common import Dir
from bioomics import BuildComplex, QueryComplex, OpenmmPdbfixer


def retrieve_pdb(args, meta):
    pdb_group = args['pdb_group']
    table_name = args['table_name']
    query = f"""
        SELECT C.pdb_id, M.chain_id, M.pdb_file
        FROM meta_chain_pdb M
        LEFT JOIN pdb_chainid C
        ON M.chain_id = C.chain_id
        LEFT JOIN meta_pdb P
        ON C.pdb_id = P.pdb_id
        WHERE P.pdb_group = {pdb_group}
            AND M.pdb_file IS NOT NULL
            AND M.chain_id NOT IN (
                SELECT chain_id
                FROM {table_name}
            )
            AND C.pdb_id IN (
                SELECT pdb_id FROM view_antibody
            )
    ;"""
    qc = QueryComplex(args['verbose'])
    rows = qc.list_data(query)
    print(f"{len(rows)} PDB files ready...")
    meta['pdb_group'] = int(pdb_group)
    meta['rows'] += len(rows)
    return rows

def run_fixer(rows, args, meta):
    for row in rows:
        pdb_id = row['pdb_id']
        pdb_endpoint = os.path.join(args['table_name'], pdb_id[:2], pdb_id)
        outdir = os.path.join(args['pdb_dir'], pdb_endpoint)
        Dir(outdir).init_dir()
        
        chain_id = row['chain_id']
        file_name = chain_id + '.pdb'
        fixed_pdb_file = os.path.join(outdir, file_name)
        if args['overwrite'] or not os.path.isfile(fixed_pdb_file):
            pdb_file = os.path.join(args['data_dir'], row['pdb_file'])
            # run pdbfixer
            op = OpenmmPdbfixer(pdb_file, args['verbose'])
            output_prefix = os.path.join(outdir, chain_id)
            info = op.fast_fixer(output_prefix)
            meta[info['status']] += 1
        else:
            meta['skip'] += 1
        if os.path.isfile(fixed_pdb_file):
            relative_path = os.path.join(pdb_endpoint, file_name)
            yield chain_id, relative_path
        else:
            meta['null'] += 1
            yield chain_id, None


if __name__ == "__main__":
    meta = defaultdict(int)
    load_dotenv()
    args = {
        'verbose': True,
        'overwrite': True,
        'pdb_dir': os.getenv('pdb_dir'),
        'data_dir': os.getenv('data_dir'),
        'chunk_size': 5,
        'table_name': 'pdb_chain_fixed',
    }
    args['pdb_group'] = sys.argv[1]
    rows = retrieve_pdb(args, meta)

    if len(rows) > 0:
        # fix pdb
        record_iter = run_fixer(rows, args, meta)
        # insertion
        bc = BuildComplex(record_iter, args['verbose'], args['chunk_size'])
        bc.insert_batch_records(args['table_name'])

    print(meta)
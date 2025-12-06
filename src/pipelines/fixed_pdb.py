#!/usr/bin/python
'''
requirements:
    - PDB raw data
    - table meta_pdb
example:
    - parallel -j1 python src/pipelines/fixed_pdb.py ::: {0..116}
    - python src/pipelines/fixed_pdb.py 0
function:
    - fix pdb files
    - insert data into table fixed_pdb
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
    # meta_pdb.pdb_file is absolute path
    query = f"""
        SELECT pdb_id, pdb_file
        FROM meta_pdb
        WHERE pdb_group = {pdb_group}
            AND pdb_id NOT IN (
                SELECT pdb_id
                FROM {table_name}
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
        pdb_endpoint = os.path.join(args['table_name'], pdb_id[:2])
        outdir = os.path.join(args['pdb_dir'], pdb_endpoint)
        output_prefix = os.path.join(outdir, pdb_id)
        fixed_pdb_file = output_prefix + '.pdb'
        if args['overwrite'] or not os.path.isfile(fixed_pdb_file):
            Dir(outdir).init_dir()
            # run pdbfixer
            op = OpenmmPdbfixer(row['pdb_file'], args['verbose'])
            info = op.fast_fixer(output_prefix)
            meta[info['status']] += 1
        else:
            meta['skip'] += 1
        if os.path.isfile(fixed_pdb_file):
            relative_path = os.path.join(pdb_endpoint, f"{pdb_id}.pdb")
            yield pdb_id, relative_path
        else:
            meta['null'] += 1
            yield pdb_id, None


if __name__ == "__main__":
    meta = defaultdict(int)
    load_dotenv()
    args = {
        'verbose': True,
        'overwrite': True,
        'pdb_dir': os.getenv('pdb_dir'),
        'chunk_size': 5,
        'table_name': 'fixed_pdb',
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
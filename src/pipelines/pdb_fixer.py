#!/usr/bin/python
'''
requirements:
    - PDB raw data
    - table meta_pdb
example:
    - parallel -j1 python src/pipelines/pdb_fixer.py ::: {0..46}
    - python src/pipelines/pdb_fixer.py 4
function:
    - fix pdb files
    - insert data into table meta_pdb_fixed
'''
from collections import defaultdict
from dotenv import load_dotenv
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import Dir, BuildComplex, QueryComplex, OpenmmPdbfixer


def retrieve_pdb(args, meta):
    pdb_group = args['pdb_group']
    # meta_pdb.pdb_file is absolute path
    query = f"""
        SELECT pdb_id, pdb_file
        FROM meta_pdb
        WHERE pdb_group = {pdb_group}
            AND pdb_id NOT IN (
                SELECT pdb_id
                FROM meta_pdb_fixed
            );
    """
    qc = QueryComplex(args['verbose'])
    rows = qc.list_data(query)
    print(f"{len(rows)} PDB files ready...")
    meta['pdb_group'] = int(pdb_group)
    meta['rows'] += len(rows)
    return rows

def run_fixer(rows, args, meta):
    for row in rows:
        pdb_id = row['pdb_id']
        outdir = os.path.join(args['pdb_dir'], 'fixed', pdb_id[:2])
        fixed_pdb_file = os.path.join(outdir, f"{pdb_id}.pdb")
        if args['overwrite'] or not os.path.isfile(fixed_pdb_file):
            Dir(outdir).init_dir()
            # run pdbfixer
            op = OpenmmPdbfixer(row['pdb_file'], args['verbose'])
            op.fast_fixer()
            status = op.save_fixer(fixed_pdb_file)
            meta[status] += 1
        else:
            meta['skip'] += 1
        if os.path.isfile(fixed_pdb_file):
            relative_path = os.path.join('fixed', pdb_id[:2], f"{pdb_id}.pdb")
            yield pdb_id, relative_path, args['pdb_group']
        else:
            meta['null'] += 1
            yield pdb_id, None, args['pdb_group']


if __name__ == "__main__":
    meta = defaultdict(int)
    load_dotenv()
    args = {
        'verbose': True,
        'overwrite': False,
        'pdb_dir': os.getenv('pdb_dir'),
        'chunk_size': 5,
    }
    args['pdb_group'] = sys.argv[1]
    rows = retrieve_pdb(args, meta)

    # fix pdb
    record_iter = run_fixer(rows, args, meta)

    # insertion
    table_name = 'meta_pdb_fixed'
    bc = BuildComplex(record_iter, args['verbose'], args['chunk_size'])
    bc.insert_meta_pdb_fixed()

    print(meta)
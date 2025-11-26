#!/usr/bin/python
'''
requirements:
    - table meta_pdb_fixed
example:
    - parallel -j2 python src/pipelines/pdb_relaxed.py ::: {0..46}
    - python src/pipelines/pdb_relaxed.py 46
function:
    - fix pdb files
    - insert data into table meta_pdb_relax
'''
from collections import defaultdict
from dotenv import load_dotenv
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import Dir, BuildComplex, QueryComplex, OpenmmSimulation


def retrieve_pdb(args, meta):
    pdb_group = args['pdb_group']
    # meta_pdb_fixed.pdb_file is relative path
    query = f"""
        SELECT pdb_id, pdb_file
        FROM meta_pdb_fixed
        WHERE pdb_group = {pdb_group}
            AND pdb_id NOT IN (
                SELECT pdb_id
                FROM meta_pdb_relaxed
            );
    """
    rows = QueryComplex(args['verbose']).list_data(query)
    print(f"{len(rows)} PDB files ready...")
    meta['pdb_group'] = int(pdb_group)
    meta['rows'] += len(rows)
    return rows

def run_simulation(rows, args, meta):
    for row in rows:
        pdb_id = row['pdb_id']
        endpoint = os.path.join('relaxed', pdb_id[:2], pdb_id)
        outdir = os.path.join(args['pdb_dir'], endpoint)
        file_name = "relaxed.pdb"
        relax_pdb_file = os.path.join(outdir, file_name)
        if args['overwrite'] or not os.path.isfile(relax_pdb_file):
            Dir(outdir).init_dir()
            # run simulation
            input_pdb = os.path.join(args['pdb_dir'], row['pdb_file'])
            status = OpenmmSimulation(input_pdb, outdir, args['verbose'])()
            meta[status] += 1
        else:
            meta['skip'] += 1
        if os.path.isfile(relax_pdb_file):
            relative_path = os.path.join(endpoint, file_name)
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
    record_iter = run_simulation(rows, args, meta)

    # insertion
    table_name = 'meta_pdb_relaxed'
    bc = BuildComplex(record_iter, args['verbose'], args['chunk_size'])
    bc.insert_meta_pdb_relaxed()

    print(meta)
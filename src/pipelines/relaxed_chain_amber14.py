#!/usr/bin/python
'''
requirements:
    - pdb fixed: table fixed_chain
example:
    - parallel -j4 python src/pipelines/relaxed_chain_amber14.py ::: {0..116}
    - python src/pipelines/relaxed_chain_amber14.py 0
function:
    - relax structure based on fixed pdb files
    - insert data into table relaxed_chain
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
    table_name = args['table_name']
    # fixed_pdb.pdb_file is relative path
    query = f"""
        SELECT P.pdb_id, F.chain_id, F.pdb_file
        FROM fixed_chain F
        LEFT JOIN meta_pdb P
        ON F.pdb_id = P.pdb_id
        WHERE P.pdb_group = {pdb_group}
            AND F.pdb_file IS NOT NULL
            AND F.chain_id NOT IN (
                SELECT chain_id FROM {table_name}
            )
            AND F.pdb_id IN (
                SELECT pdb_id from view_antibody
            );
    """
    rows = QueryComplex(args['verbose']).list_data(query)
    print(f"{len(rows)} PDB files ready...")
    meta['pdb_group'] = int(pdb_group)
    meta['rows'] += len(rows)
    return rows

def run_simulation(rows, args, meta):
    for row in rows:
        pdb_id, chain_id = row['pdb_id'], row['chain_id'], 
        original_energy, relaxed_energy = None, None
        endpoint = os.path.join(args['table_name'], pdb_id[:2], pdb_id, chain_id)
        outdir = os.path.join(args['pdb_dir'], endpoint)
        Dir(outdir).init_dir()
        
        # run simulation
        input_pdb = os.path.join(args['pdb_dir'], row['pdb_file'])
        curr_args = {
            'pdb_id': pdb_id,
            'outdir': outdir,
        }
        curr_args.update(args)
        info = OpenmmSimulation(input_pdb, curr_args).relax_pdb()
        meta[info['status']] += 1
        original_energy = info.get('original_energy')
        relaxed_energy = info.get('relaxed_energy')
        
        # build record
        rec0 = (row['chain_id'], original_energy, relaxed_energy,)
        relax_pdb_file = os.path.join(outdir, 'relaxed.pdb')
        if os.path.isfile(relax_pdb_file):
            yield rec0 + (endpoint,)
        else:
            meta['null'] += 1
            yield rec0 + (None,)

if __name__ == "__main__":
    meta = defaultdict(int)
    load_dotenv()
    args = {
        'verbose': True,
        'pdb_dir': os.getenv('pdb_dir'),
        'chunk_size': 5,
        'table_name': 'relaxed_chain_amber14',
        'force_fields': ['amber14-all.xml', 'amber14/tip3pfb.xml'],
        'max_iterations': 0,
        'steps': 100_000,
    }
    args['pdb_group'] = sys.argv[1]
    rows = retrieve_pdb(args, meta)

    if len(rows) > 0:
        # relax pdb
        record_iter = run_simulation(rows, args, meta)
        # insertion
        bc = BuildComplex(record_iter, args['verbose'], args['chunk_size'])
        bc.insert_batch_records(args['table_name'])

    print(meta)
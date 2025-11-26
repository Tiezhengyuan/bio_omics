#!/usr/bin/python
'''
requirements:
    - table meta_pdb
example:
    - parallel -j2 python src/pipelines/pdb_chain_type.py ::: {0..46}
    - python src/pipelines/pdb_chain_type.py 46
function:
    - fix pdb files
    - insert data into table chain_type
'''
from collections import defaultdict
from dotenv import load_dotenv
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import Pdb, BuildComplex, QueryComplex


def retrieve_pdb(args, meta):
    pdb_group = args['pdb_group']
    query = f"""
        SELECT pdb_id, pdb_file
        FROM meta_pdb
        WHERE pdb_group = {pdb_group}
            AND pdb_id NOT IN (
                SELECT pdb_id
                FROM pdb_chain_type
            );
    """
    rows = QueryComplex(args['verbose']).list_data(query)
    print(f"{len(rows)} PDB files ready...")
    meta['pdb_group'] = int(pdb_group)
    meta['rows'] += len(rows)
    return rows

def run(rows, args, meta):
    for row in rows:
        meta['pdb_record'] += 1
        p = Pdb(row['pdb_file'], None, args['verbose'])
        p.load_structure()
        chains = p.parse_chains()
        if chains:
            for chain_id, pdb_id, model_no, chain_no, chain in chains:
                counts = p.residue_type(chain)
                rec = (pdb_id, chain_id, counts['type'], counts['total_res'], \
                    counts['non_water_res'],  counts['aa_res'], counts['dna_res'], \
                    counts['rna_res'], counts['water_res'], counts['ion_res'], \
                    counts['hetero_res'], counts['other_res'])
                meta['chain'] += 1
                yield rec
        else:
            rec0 = (row["pdb_id"],) + (None,) * 11
            yield rec0

if __name__ == "__main__":
    meta = defaultdict(int)
    load_dotenv()
    args = {
        'verbose': True,
        'chunk_size': 2,
    }
    args['pdb_group'] = sys.argv[1]
    rows = retrieve_pdb(args, meta)

    # fix pdb
    record_iter = run(rows, args, meta)

    # insertion
    table_name = 'meta_pdb_relaxed'
    bc = BuildComplex(record_iter, args['verbose'], args['chunk_size'])
    bc.insert_pdb_chain_type()

    print(meta)
#!/usr/bin/python
'''
requirements:
    - table meta_pdb, pdb_chainid
example:
    - parallel -j4 python app.py pdb_chain_type ::: {0..237}
    - python app.py pdb_chain_type 0
function:
    - fix pdb files
    - insert data into table pdb_chain_type
'''

from src.helper import *

@header
def retrieve_pdb(params, meta):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        SELECT pdb_id, relative_path
        FROM meta_pdb
        WHERE pdb_group = {pdb_group}
            AND pdb_id NOT IN (
                SELECT pdb_id FROM {table_name}
            );
    """
    rows = QueryComplex(params['verbose']).list_data(query)
    print(f"{len(rows)} PDB files ready...")
    meta['pdb_group'] = int(pdb_group)
    meta['rows'] += len(rows)
    for row in rows:
        pdb_file = os.path.join(params['pdb_dir'], row['relative_path'])
        if os.path.isfile(pdb_file):
            row['pdb_file'] = pdb_file
            yield row

def detect_chain_type(rows, params, meta):
    for row in rows:
        pdb_id = row['pdb_id']
        meta['pdb_record'] += 1
        p = Pdb(row['pdb_file'], None, params['verbose'])
        p.load_structure()

        chains = p.parse_chains()
        if chains:
            for item in chains:
                counts = p.residue_type(item['chain'])
                rec = (pdb_id, item['chain_id'], counts['type'], counts['total_res'], \
                    counts['non_water_res'],  counts['aa_res'], counts['dna_res'], \
                    counts['rna_res'], counts['water_res'], counts['ion_res'], \
                    counts['hetero_res'], counts['other_res'])
                meta['chain'] += 1
                yield rec
        else:
            rec0 = (pdb_id,) + (None,) * 11
            yield rec0


if __name__ == "__main__":
    params.update({
        'chunk_size': 50,
        'table_name': 'pdb_chain_type',
    })
    rows = retrieve_pdb(params, meta)

    # retrieve chain types
    record_iter = detect_chain_type(rows, params, meta)
    bc = BuildComplex(record_iter, params['verbose'], params['chunk_size'])
    bc.insert_batch_records(params['table_name'])

    footer(meta)
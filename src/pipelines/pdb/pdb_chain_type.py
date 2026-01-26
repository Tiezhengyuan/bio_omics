#!/usr/bin/python
'''
requirements:
    - table meta_pdb, pdb_chainid
example:
    - parallel -j4 python src/pipelines/pdb_chain_type.py ::: {0..117}
    - python src/pipelines/pdb_chain_type.py 0
function:
    - fix pdb files
    - insert data into table pdb_chain_type
'''

from helper import *

def retrieve_pdb(args, meta):
    pdb_group = args['pdb_group']
    table_name = args['table_name']
    query = f"""
        SELECT pdb_id, relative_path
        FROM meta_pdb
        WHERE pdb_group = {pdb_group}
            AND pdb_id NOT IN (
                SELECT pdb_id FROM {table_name}
            );
    """
    rows = QueryComplex(args['verbose']).list_data(query)
    print(f"{len(rows)} PDB files ready...")
    meta['pdb_group'] = int(pdb_group)
    meta['rows'] += len(rows)
    for row in rows:
        pdb_file = os.path.join(args['pdb_dir'], row['relative_path'])
        if os.path.isfile(pdb_file):
            row['pdb_file'] = pdb_file
            yield row

def detect_chain_type(rows, args, meta):
    for row in rows:
        pdb_id = row['pdb_id']
        meta['pdb_record'] += 1
        p = Pdb(row['pdb_file'], None, args['verbose'])
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
    args = {
        'verbose': True,
        'chunk_size': 50,
        'pdb_dir': os.getenv('pdb_dir'),
        'table_name': 'pdb_chain_type',
    }
    args['pdb_group'] = sys.argv[1]
    rows = retrieve_pdb(args, meta)

    # retrieve chain types
    record_iter = detect_chain_type(rows, args, meta)
    bc = BuildComplex(record_iter, args['verbose'], args['chunk_size'])
    bc.insert_batch_records(args['table_name'])

    print(meta)
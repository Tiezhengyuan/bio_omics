#!/usr/bin/python
'''
requirements:
    - meta_pdb
    - fixed_pdb
example:
    - parallel -j4 python src/pipelines/fixed_chainid.py ::: {0..116}
    - python src/pipelines/fixed_chainid.py 0
function:
    - parse chain.id between original and fixed PDB
'''

from helper import *

def retrieve_data(pdb_group, args, meta):
    table_name = args['table_name']
    query = f"""
        SELECT M.pdb_id, M.relative_path, F.relative_path AS fixed_relative
        FROM meta_pdb M, fixed_pdb F
        WHERE M.pdb_id = F.pdb_id
            AND M.pdb_group = {pdb_group}
            AND M.relative_path IS NOT NULL
            AND F.relative_path IS NOT NULL
            AND M.pdb_id NOT IN (
                SELECT pdb_id FROM {table_name}
            )
    ;"""
    rows = QueryComplex(args['verbose']).list_data(query)
    meta['pdb_group'] = int(pdb_group)
    meta['rows'] = len(rows)
    return rows

def parse(rows, args, meta):
    for row in rows:
        pdb_id = row['pdb_id']
        try:
            # original chains
            input_pdb = os.path.join(args['pdb_dir'], row['relative_path'])
            p1 = Pdb(input_pdb, None, args['verbose'])
            p1.load_structure()
            # fixed chains
            fixed_pdb = os.path.join(args['pdb_dir'], row['fixed_relative'])
            p2 = Pdb(fixed_pdb, None, args['verbose'])
            p2.load_structure()
            
            # parse
            for p1_model, p2_model in zip(p1.structure, p2.structure):
                p1_model_no, p2_model_no = p1_model.id, p2_model.id
                for p1_chain, p2_chain in zip(p1_model, p2_model):
                    p1_chain_no, p2_chain_no = p1_chain.id, p2_chain.id
                    meta['num_record'] += 1
                    rec = (pdb_id, p1_model_no, p2_model_no, p1_chain_no, p2_chain_no)
                    if args['verbose']:
                        print('Insertion record: ', rec)
                    yield rec
        except Exception as e:
            print(e)
            meta['fail'] += 1

##################################
if __name__ == "__main__":
    pdb_group = sys.argv[1]
    args = {
        'verbose': True,
        'overwrite': True,
        'pdb_dir': os.getenv('pdb_dir'),
        'data_dir': os.getenv('data_dir'),
        'chunk_size': 50,
        'table_name': 'fixed_chainid',
    }
    rows = retrieve_data(pdb_group, args, meta)

    if len(rows) > 0:
        print(f"Try to insert {len(rows)} data into {args['table_name']}...")
        record_iter = parse(rows, args, meta)
        bc = BuildComplex(record_iter, args['verbose'], args['chunk_size'])
        bc.insert_batch_records(args['table_name'])

    print(meta)


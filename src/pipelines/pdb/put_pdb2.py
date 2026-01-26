#!/usr/bin/python
'''
requirements:
    - PDB raw data should be downloaded to the local
example:
    - python src/pipelines/put_pdb.py
function:
    - build table pdb and insert data into table pdb
'''
from dotenv import load_dotenv
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import PullPdb, BuildComplex, QueryComplex


if __name__ == "__main__":
    load_dotenv()
    args = {
        'verbose': True,
        'new_pdb': False,
        'pdb_dir': os.getenv('pdb_dir'),
        'data_format': 'mmCIF',
        'chunk_size': 100,
    }

    # table pdb is foundamental table
    pdb_dir = os.path.join(args['pdb_dir'], args['data_format'])
    table_name = 'pdb'
    if args['new_pdb']:
        pdb_path_dict = PullPdb(pdb_dir).get_raw_pdb()
        pdb_files = pdb_path_dict.values()
        bc = BuildComplex(pdb_files, args['verbose'], args['chunk_size'])
        bc.empty_table(table_name)
        bc.insert_pdb()
    else:
        db_pdb_ids = QueryComplex(args['verbose']).unique_values(table_name, 'pdb_id')
        pdb_path_dict = PullPdb(pdb_dir).filter_out_raw_pdb(db_pdb_ids)
        pdb_files = pdb_path_dict.values()
        if pdb_files:
            print(f"Try to insert {len(pdb_files)} PDB data into table {table_name}...")
            bc = BuildComplex(pdb_files, args['verbose'], args['chunk_size'])
            bc.insert_pdb()
        else:
            print(f'Synchronization of PDB with {table_name} is good.')

    
    print("Done. Great\n\n")

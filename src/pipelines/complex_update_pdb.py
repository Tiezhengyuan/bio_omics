'''
example: python src/pipelines/complex_pdb_update.py
synchronize PDB/UniProt/alignments with data in database COMPLEX.
'''
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import QueryComplex,  BuildComplex, PullPdb

##################################
def update_pdb(args):
    qc = QueryComplex(args['verbose'])
    table_name = 'pdb'
    db_pdb_ids = qc.get_unique_values(table_name, 'pdb_id')
    pdb_path_dict = PullPdb(args['pdb_dir']).filter_out_raw_pdb(db_pdb_ids)
    pdb_files = pdb_path_dict.values()
    if pdb_files:
        print(f"Try to insert {len(pdb_files)} PDB data into table {table_name}...")
        BuildComplex(pdb_files, args['verbose']).insert_pdb()
    else:
        print(f'Synchronization of PDB with {table_name} is good.')

def update_compound(args):
    qc = QueryComplex(args['verbose'])
    table_name = 'compound'
    db_pdb_ids = qc.get_unique_values(table_name, 'pdb_id')
    pdb_path_dict = PullPdb(args['pdb_dir']).filter_out_raw_pdb(db_pdb_ids)
    pdb_files = pdb_path_dict.values()
    if pdb_files:
        print(f"Try to insert {len(pdb_files)} PDB data into table {table_name}...")
        BuildComplex(pdb_files, args['verbose']).insert_compound()
    else:
        print(f'Synchronization of PDB with {table_name} is good.')

def update_chain(args):
    qc = QueryComplex(args['verbose'])
    table_name = 'chain'
    db_pdb_ids = qc.get_unique_values(table_name, 'pdb_id')
    pdb_path_dict = PullPdb(args['pdb_dir']).filter_out_raw_pdb(db_pdb_ids)
    pdb_files = pdb_path_dict.values()
    if pdb_files:
        print(f"Try to insert {len(pdb_files)} PDB data into table {table_name}...")
        BuildComplex(pdb_files, args['verbose']).insert_chain(check=True)
    else:
        print(f'Synchronization of PDB with {table_name} is good.')

def update_keyword(args):
    qc = QueryComplex(args['verbose'])
    table_name = 'keyword'
    db_pdb_ids = qc.get_unique_values(table_name, 'pdb_id')
    pdb_path_dict = PullPdb(args['pdb_dir']).filter_out_raw_pdb(db_pdb_ids)
    pdb_files = pdb_path_dict.values()
    if pdb_files:
        print(f"Try to insert {len(pdb_files)} PDB data into table {table_name}...")
        BuildComplex(pdb_files, args['verbose']).insert_keyword()
    else:
        print(f'Synchronization of PDB with {table_name} is good.')

##################################
if __name__ == "__main__":
    args = {
        'verbose': True,
        'pdb_dir': '/home/yuan/data/pdb/rsync',
        'chunk_size': 20,
    }
    # note: new PDB data should be synchronized to the local.
    # pdb related tables
    update_pdb(args)
    update_compound(args)
    update_chain(args)
    update_keyword(args)

'''
note:
    - PDB files should be synchronized to the local.
    - The table pdb is ready.
example: python src/pipelines/put_chain.py
function:
    synchronize PDB/UniProt/alignments with data in database COMPLEX.
'''
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import QueryComplex,  BuildComplex

##################################
if __name__ == "__main__":
    args = {
        'verbose': True,
        'pdb_dir': '/home/yuan/data/pdb/rsync',
        'chunk_size': 500,
    }

    qc = QueryComplex(args['verbose'])
    table_name = 'pdb_chain'
    pdb_files = qc.newly_pdb_ids(table_name)
    if pdb_files:
        print(f"Try to insert {len(pdb_files)} PDB data into table {table_name}...")
        bc = BuildComplex(pdb_files, args['verbose'], args['chunk_size'])
        bc.insert_pdb_chain()
    else:
        print(f'Synchronization of PDB with {table_name} is good.')


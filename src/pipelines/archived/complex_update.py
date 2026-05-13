'''
example: python src/pipelines/complex_update.py
synchronize PDB/UniProt/alignments with data in database COMPLEX.
Note: those scripts should be executed precdently.
    - pdb_blast.py
    - pdb_igblast.py
    - complex_pdb_insert.py
    - complex_uniprot_insert.py
'''
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import BuildComplex, FilterDb

##################################

def update_align_uniprot(args):
    table_name = 'align_uniprot'
    print(f"Try to insert UniProt data into table {table_name}...")
    fd = FilterDb(args['data_dir'], args['verbose'])
    aln_iter = fd.uniprot_align_check_iter(table_name)
    bc = BuildComplex(aln_iter, args['verbose'], args['chunk_size'])
    bc.insert_align_uniprot()
    print(f'Synchronization of PDB with {table_name} is good.')

def update_align_imgt(args):
    table_name = 'align_imgt'
    print(f"Try to insert UniProt data into table {table_name}...")
    fd = FilterDb(args['data_dir'], args['verbose'])
    aln_iter = fd.igblastp_align_check_iter(table_name)
    bc = BuildComplex(aln_iter, args['verbose'], args['chunk_size'])
    bc.insert_align_imgt()
    print(f'Synchronization of PDB with {table_name} is good.')

def update_align_vregion(args):
    table_name = 'align_vregion'
    print(f"Try to insert UniProt data into table {table_name}...")
    fd = FilterDb(args['data_dir'], args['verbose'])
    aln_iter = fd.igblastp_align_check_iter(table_name)
    bc = BuildComplex(aln_iter, args['verbose'])
    bc.insert_align_vregion()
    print(f'Synchronization of PDB with {table_name} is good.')

##################################
if __name__ == "__main__":
    args = {
        'verbose': True,
        'pdb_dir': '/home/yuan/data/pdb/rsync',
        'data_dir': '/home/yuan/data/pdb/data',
        'chunk_size': 20,
    }
    # update blastp results
    update_align_uniprot(args)
    # update igblastp results
    update_align_imgt(args)
    update_align_vregion(args)


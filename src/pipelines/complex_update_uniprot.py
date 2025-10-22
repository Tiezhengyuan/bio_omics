'''
example: python src/pipelines/db_complex_update.py
synchronize PDB/UniProt/alignments with data in database COMPLEX.
'''
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import QueryComplex,  BuildComplex, Uniprot


def update_uniprot(args):
    qc = QueryComplex(args['verbose'])
    table_name = 'uniprot'
    db_ids = qc.get_unique_values(table_name, 'uniprot_acc')
    filtered_records_iter = Uniprot(args['uniprot_sprot_dat']).filter_dat_iter(db_ids)
    print(f"Try to insert UniProt data into table {table_name}...")
    bc = BuildComplex(filtered_records_iter, args['verbose'])
    bc.insert_uniprot()
    print(f'Synchronization of PDB with {table_name} is good.')

def update_uniprot_keyword(args):
    qc = QueryComplex(args['verbose'])
    table_name = 'uniprot_keyword'
    db_ids = qc.get_unique_values(table_name, 'uniprot_acc')
    filtered_records_iter = Uniprot(args['uniprot_sprot_dat']).filter_dat_iter(db_ids)
    print(f"Try to insert UniProt data into table {table_name}...")
    bc = BuildComplex(filtered_records_iter, args['verbose'])
    bc.insert_uniprot_keyword()
    print(f'Synchronization of PDB with {table_name} is good.')

def update_uniprot_dbxrefs(args):
    qc = QueryComplex(args['verbose'])
    table_name = 'uniprot_dbxrefs'
    db_ids = qc.get_unique_values(table_name, 'uniprot_acc')
    filtered_records_iter = Uniprot(args['uniprot_sprot_dat']).filter_dat_iter(db_ids)
    print(f"Try to insert UniProt data into table {table_name}...")
    bc = BuildComplex(filtered_records_iter, args['verbose'])
    bc.insert_uniprot_dbxrefs()
    print(f'Synchronization of PDB with {table_name} is good.')

##################################
if __name__ == "__main__":
    args = {
        'verbose': True,
        'pdb_dir': '/home/yuan/data/pdb/rsync',
        'uniprot_sprot_dat': '/home/yuan/data/UniProt/uniprot_sprot.dat.gz',
        'data_dir': '/home/yuan/data/pdb/data',
    }

    # needed when new uniprot data is available
    # uniprot related tables
    update_uniprot(args)
    update_uniprot_keyword(args)
    update_uniprot_dbxrefs(args)


'''
example: python src/pipelines/complex_update.py
synchronize PDB/UniProt/alignments with data in database COMPLEX.
'''
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import QueryComplex,  BuildComplex, PullPdb, FilterDb

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
        'uniprot_sprot_dat': '/home/yuan/data/UniProt/uniprot_sprot.dat.gz',
        'data_dir': '/home/yuan/data/pdb/data',
        'chunk_size': 20,
    }
    # needed when new PDB data available
    # pdb related tables
    # update_pdb(args)
    # update_compound(args)
    # update_chain(args)
    # update_keyword(args)

    # needed when pdb_analyze.py is executed
    # update_align_uniprot(args)
    # update_align_imgt(args)
    update_align_vregion(args)


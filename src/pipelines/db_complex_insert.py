#!/usr/bin/python
'''
example: python src/pipelines/db_complex_insert.py
function: build table pdb and insert data into table pdb
'''
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)


from bioomics import PullPdb, BuildComplex, QueryComplex, \
    QueryDb, FilterDb, Uniprot, Dir


def build_pdb_tables(args):
    pdb_path_dict = PullPdb(args['pdb_dir']).get_raw_pdb()
    pdb_files = pdb_path_dict.values()
    bc = BuildComplex(pdb_path_files, False)
    # table pdb
    bc.insert_pdb()
    # table compound
    bc.insert_compound()
    # table chain
    bc.insert_chain()
    # table pdb_keyword
    bc.insert_pdb_keyword()

def build_uniprot_tables(args):
    '''
    Data comes from one dat file.
    Don't recommend updates. Delete all rows before insertion 
    '''
    # table uniprot
    # records_iter = Uniprot(args['uniprot_sprot_dat']).dat_iter()
    # BuildComplex(records_iter, True).insert_uniprot()

    # table uniprot_keyword
    # records_iter = Uniprot(args['uniprot_sprot_dat']).dat_iter()
    # BuildComplex(records_iter, True).insert_uniprot_keyword()

    # table uniprot_dbxrefs
    # records_iter = Uniprot(args['uniprot_sprot_dat']).dat_iter()
    # BuildComplex(records_iter, True).insert_uniprot_dbxrefs()

def build_blastp_table(args):
    # table uniprot_aling
    # collect aln files of which data not in DB
    aln_iter = FilterDb(args['data_dir']).uniprot_align()
    BuildComplex(aln_iter, True).insert_uniprot_align()


if __name__ == "__main__":
    args = {
        'pdb_dir': '/home/yuan/data/pdb/rsync',
        'uniprot_sprot_dat': '/home/yuan/data/UniProt/uniprot_sprot.dat.gz',
        'data_dir': '/home/yuan/data/pdb/data',
    }
    
    # scan pdb files
    # build_pdb_tables(args)
    
    # data are from *.dat
    # build_uniprot_tables(args)

    # 
    build_blastp_table(args)
    
    print("Done. Great\n\n")

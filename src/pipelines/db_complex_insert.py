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
    bc = BuildComplex(pdb_files, args['verbose'])
    # table pdb
    bc.insert_pdb()
    # table compound
    bc.insert_compound()
    # table chain
    bc.insert_chain()
    # table keyword
    bc.insert_keyword()

def build_uniprot_tables(args):
    '''
    Data comes from one dat file.
    Don't recommend updates. Delete all rows before insertion 
    '''
    # table uniprot
    records_iter = Uniprot(args['uniprot_sprot_dat']).dat_iter()
    BuildComplex(records_iter, args['verbose']).insert_uniprot()

    # table uniprot_keyword
    records_iter = Uniprot(args['uniprot_sprot_dat']).dat_iter()
    BuildComplex(records_iter, args['verbose']).insert_uniprot_keyword()

    # table uniprot_dbxrefs
    records_iter = Uniprot(args['uniprot_sprot_dat']).dat_iter()
    BuildComplex(records_iter, args['verbose']).insert_uniprot_dbxrefs()

def build_blast(args):
    # table align_uniprot
    # collect aln files of which data not in DB
    aln_iter = FilterDb(args['data_dir']).uniprot_align_iter()
    bc = BuildComplex(aln_iter, args['verbose'])
    bc.empty_table('align_uniprot')
    bc.insert_align_uniprot()

def build_igblast(args):
    # table align_imgt
    aln_fa_iter = FilterDb(args['data_dir']).igblastp_align_iter()
    bc = BuildComplex(aln_fa_iter, args['verbose'])
    bc.empty_table('align_imgt')
    bc.insert_align_imgt()

    # table align_vregion
    aln_fa_iter = FilterDb(args['data_dir']).igblastp_align_iter()
    bc = BuildComplex(aln_fa_iter, args['verbose'])
    bc.empty_table('align_vregion')
    bc.insert_align_vregion()

if __name__ == "__main__":
    args = {
        'verbose': True,
        'pdb_dir': '/home/yuan/data/pdb/rsync',
        'uniprot_sprot_dat': '/home/yuan/data/UniProt/uniprot_sprot.dat.gz',
        'data_dir': '/home/yuan/data/pdb/data',
    }
    
    # scan pdb files
    # build_pdb_tables(args)

    # alignments
    # build_align_uniprot(args)

    # data are from *.dat
    # build_blast(args)
    build_igblast(args)
    
    print("Done. Great\n\n")

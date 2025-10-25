#!/usr/bin/python
'''
example: python src/pipelines/complex_insert.py
function: build table pdb and insert data into table pdb
'''
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import BuildComplex, Uniprot

def insert(args):
    '''
    Data comes from one dat file.
    Don't recommend updates. Delete all rows before insertion 
    '''
    # table uniprot
    records_iter = Uniprot(args['uniprot_sprot_dat']).dat_iter()
    bc = BuildComplex(records_iter, args['verbose'])
    bc.empty_table('uniprot')
    bc.insert_uniprot()

    # table uniprot_keyword
    records_iter = Uniprot(args['uniprot_sprot_dat']).dat_iter()
    bc = BuildComplex(records_iter, args['verbose'])
    bc.empty_table('uniprot_keyword')
    bc.insert_uniprot_keyword()

    # table uniprot_dbxrefs
    records_iter = Uniprot(args['uniprot_sprot_dat']).dat_iter()
    bc = BuildComplex(records_iter, args['verbose'])
    bc.empty_table('uniprot_dbxrefs')
    bc.insert_uniprot_dbxrefs()

if __name__ == "__main__":
    args = {
        'verbose': True,
        'uniprot_sprot_dat': '/home/yuan/data/UniProt/uniprot_sprot.dat.gz',
    }
    print('Insert data into tables...')
    insert(args)

    print("Done. Great\n\n")

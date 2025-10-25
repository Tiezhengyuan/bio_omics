#!/usr/bin/python
'''
example: python src/pipelines/complex_pdb_insert.py
function: build table pdb and insert data into table pdb
'''
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import PullPdb, BuildComplex

def insert(args):
    pdb_path_dict = PullPdb(args['pdb_dir']).get_raw_pdb()
    pdb_files = pdb_path_dict.values()
    bc = BuildComplex(pdb_files, args['verbose'])
    # table pdb is foundamental table
    bc.empty_table('pdb')
    bc.insert_pdb()
    # table compound
    bc.empty_table('compound')
    bc.insert_compound()
    # table chain
    bc.empty_table('chain')
    bc.insert_chain()
    # table keyword
    bc.empty_table('keyword')
    bc.insert_keyword()


if __name__ == "__main__":
    args = {
        'verbose': True,
        'pdb_dir': '/home/yuan/data/pdb/rsync',
    }
    # insert data into empty tables
    insert(args)

    
    print("Done. Great\n\n")

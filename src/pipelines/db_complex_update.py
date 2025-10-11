'''
build table pdb and insert data into table pdb
'''
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import PullPdb, BuildComplex, QueryDb, Uniprot

##################################
def update_pdb_tables(args):
    pass

def update_uniprot_tables(args):
    pass

##################################
if __name__ == "__main__":
    args = {
        'pdb_dir': '/home/yuan/data/pdb/rsync',
    }
    # pdb related tables
    update_pdb_tables(args)

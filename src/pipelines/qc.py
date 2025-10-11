#!/usr/bin/python
'''
example: python src/pipelines/qc.py
functions: data quality control
'''
from Bio import SeqIO
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)


from bioomics import Dir, PullPdb


def dup_raw_pdb(args):
    print(f"\nDetect duplicated PDB in {args['pdb_dir']}...")
    res = {}
    pdb_iter = PullPdb(args['pdb_dir']).iter_raw_pdb()
    for pdb_id, pdb_file in pdb_iter:
        if pdb_id not in res:
            res[pdb_id] = pdb_file
        else:
            print(pdb_file, res[pdb_id])

def dup_data_pdb(args):
    print(f"\nDetect duplicated PDB label in {args['data_dir']}...")
    res = {}
    dir_iter = Dir(args['data_dir']).recursive_dirs()
    for path in dir_iter:
        basename = os.path.basename(path)
        if basename not in res:
            res[basename] = path
        else:
            print(path, res[basename])


if __name__ == "__main__":
    args = {
        'pdb_dir': '/home/yuan/data/pdb/rsync',
        'data_dir': '/home/yuan/data/pdb/data',
    }

    dup_raw_pdb(args)
    dup_data_pdb(args)
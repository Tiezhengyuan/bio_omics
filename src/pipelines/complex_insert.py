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

from bioomics import Dir, BuildComplex, FilterDb


# table align_uniprot
def build_blast(args):
    # collect aln files of which data not in DB
    aln_iter = FilterDb(args['data_dir']).uniprot_align_iter()
    bc = BuildComplex(aln_iter, args['verbose'])
    bc.empty_table('align_uniprot')
    bc.insert_align_uniprot()

# table align_imgt
def build_igblast(args):
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
        'data_dir': '/home/yuan/data/pdb/data',
        'alphafold_model': '/home/yuan/rawdata/alphafold.ebi.ac.uk',
    }
    # build_blast(args)
    # build_igblast(args)
    

    print("Done. Great\n\n")


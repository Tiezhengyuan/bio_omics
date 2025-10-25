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

def build_modelling(args):
    records = []
    indir = os.path.join(args['alphafold_model'], 'swissprot_pdb_v6')
    file_iter = Dir(indir).recursive_files()
    for path in file_iter:
        file_name = os.path.basename(path)
        uniprot_acc = file_name.split('-')[1]
        rec = {
            'uniprot_acc': uniprot_acc,
            'model': 'alphafold',
            'pdb_path': path,
        }
        records.append(rec)

    bc = BuildComplex(records, args['verbose'])
    table_name = 'uniprot_structure'
    bc.empty_table(table_name)
    bc.insert_uniprot_structure()

if __name__ == "__main__":
    args = {
        'verbose': True,
        'data_dir': '/home/yuan/data/pdb/data',
        'alphafold_model': '/home/yuan/rawdata/alphafold.ebi.ac.uk',
    }
    # build_blast(args)
    # build_igblast(args)
    
    # alphafold database should be downloaded
    build_modelling(args)

    print("Done. Great\n\n")

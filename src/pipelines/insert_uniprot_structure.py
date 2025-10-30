#!/usr/bin/python
'''
requirements:
    - alphafold database should be downloaded
    - table uniprot
example:
    - python src/pipelines/insert_uniprot_structure.py
function:
    - build table pdb and insert data into table pdb
'''
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import Dir, BuildComplex


if __name__ == "__main__":
    args = {
        'verbose': True,
        'data_dir': '/home/yuan/data/pdb/data',
        'alphafold_model': '/home/yuan/rawdata/alphafold.ebi.ac.uk',
    }

    # collect data
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

    # insertion
    bc = BuildComplex(records, args['verbose'])
    table_name = 'uniprot_structure'
    bc.empty_table(table_name)
    bc.insert_uniprot_structure()

    print("Done. Great\n\n")


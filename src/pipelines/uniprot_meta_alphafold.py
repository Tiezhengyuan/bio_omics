#!/usr/bin/python
'''
requirements:
    - alphafold database should be downloaded
    - table uniprot
example:
    - python src/pipelines/uniprot_meta_alphafold.py
function:
    - put data into table uniprot_meta_alphafold
'''
from collections import defaultdict
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import Dir, BuildComplex, QueryComplex

def search_data(meta):
    query = """
        SELECT uniprot_acc
        FROM uniprot
        WHERE uniprot_acc NOT IN (
            SELECT uniprot_acc
            FROM uniprot_meta_alphafold
        );
    """
    rows = QueryComplex(args['verbose']).list_data(query)
    meta['rows'] = len(rows)
    ids = [i['uniprot_acc'] for i in rows]
    return ids

def build_iter(file_iter, ids, meta):
    for path in file_iter:
        file_name = os.path.basename(path)
        uniprot_acc = file_name.split('-')[1]
        if uniprot_acc in ids:
            rec = (uniprot_acc, path,)
            meta['records'] += 1
            yield rec

if __name__ == "__main__":
    meta = defaultdict(int)
    args = {
        'verbose': True,
        'data_dir': '/home/yuan/data/pdb/data',
        'alphafold_model': '/home/yuan/rawdata/alphafold.ebi.ac.uk/swissprot_pdb_v6',
    }
    # collect local path
    file_iter = Dir(args['alphafold_model']).recursive_files()
    # collect uniprot_acc not in table uniprot_meta_alphafold
    table_name = 'uniprot_structure'
    ids = search_data(meta)
    # collect data
    record_iter = build_iter(file_iter, ids, meta)

    # insertion
    print(f"Try to insert newly records into {table_name}...")
    bc = BuildComplex(record_iter, args['verbose'])
    bc.insert_uniprot_meta_alphafold()

    print(meta)
    print("Done. Great\n\n")


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

from helper import *

def search_data(args, meta):
    table_name = args['table_name']
    query = f"""
        SELECT uniprot_acc
        FROM uniprot
        WHERE uniprot_acc NOT IN (
            SELECT uniprot_acc
            FROM {table_name}
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
    args = {
        'verbose': True,
        'data_dir': os.getenv('data_dir'),
        'alphafold_model': '/home/yuan/rawdata/alphafold.ebi.ac.uk/swissprot_pdb_v6',
        'table_name': 'uniprot_meta_alphafold',
    }
    # collect local path
    file_iter = Dir(args['alphafold_model']).recursive_files()
    # collect uniprot_acc not in table uniprot_meta_alphafold
    table_name = 'uniprot_structure'
    ids = search_data(args, meta)
    # collect data
    record_iter = build_iter(file_iter, ids, meta)

    # insertion
    print(f"Try to insert newly records into {table_name}...")
    bc = BuildComplex(record_iter, args['verbose'])
    bc.insert_batch_records(args['table_name'])

    print(meta)
    print("Done. Great\n\n")


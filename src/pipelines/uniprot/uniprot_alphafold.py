#!/usr/bin/python
'''
requirements:
    - alphafold database should be downloaded
    - table uniprot
example:
    - python app.py uniprot_alphafold
function:
    - put path of predict pdb into uniprot_alphafold
'''

from src.helper import *
from bioomics import Dir

def search_data(params, meta):
    table_name = params['table_name']
    query = f"""
        SELECT uniprot_acc FROM uniprot
        WHERE uniprot_acc NOT IN (
            SELECT uniprot_acc FROM {table_name}
        );
    """
    rows = QueryComplex(params['verbose']).list_data(query)
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
    params.update({
        'alphafold_model': '/home/yuan/rawdata/alphafold.ebi.ac.uk/swissprot_pdb_v6',
        'table_name': 'uniprot_alphafold',
    })

    # collect local path
    file_iter = Dir(params['alphafold_model']).recursive_files()
    # collect uniprot_acc not in the table
    table_name = 'uniprot_structure'
    ids = search_data(params, meta)
    # collect data
    record_iter = build_iter(file_iter, ids, meta)

    print(f"Try to insert newly records into {table_name}...")
    bc = BuildComplex(params['verbose'])
    bc.insert_batch_records(record_iter, params['table_name'])

    footer(meta)

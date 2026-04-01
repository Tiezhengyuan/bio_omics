#!/usr/bin/python
'''
requirements:
    - meta_combo2, meta_combo3
    - table view_antibody
example:
    - python app.py pdb_complex
function:
    - identify complex
    - insert data into table pdb_complex
'''
from src.helper import *


def antibody_antigen(args, meta):
    res = QueryComplex(params['verbose']).abag_combo2()
    for v in res.values():
        rec = (v['pdb_id'], v['combo_id'], v['chain_combo'], 
            v['chain_label'], 2, 'AbAg')
        yield rec

def heavy_light_antigen(args, meta):
    res = QueryComplex(params['verbose']).abag_combo3()
    for v in res.values():
        rec = (v['pdb_id'], v['combo_id'], v['chain_combo'], 
            v['chain_label'], 3, 'HLAg')
        yield rec


if __name__ == "__main__":
    params = {
        'verbose': True,
        'table_name': 'pdb_complex',
    }
    UpdateComplex(params['table_name'], params['verbose']).empty_tables()

    # antibody complex
    pool = [antibody_antigen, heavy_light_antigen]
    for func in pool:
        record_iter = func(params, meta)
        bc = BuildComplex(record_iter, params['verbose'])
        bc.insert_batch_records(params['table_name'])

    footer(params, meta)
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
from bioomics import QueryAntibody

# dimer
def heav_light(meta):
    res = QueryAntibody(params['verbose']).combo2_heavy_light()
    for v in res.values():
        rec = (v['pdb_id'], v['combo_id'], v['chain_combo'], 
            v['chain_label'], 2, 'antibody')
        meta['heavy-light'] += 1
        yield rec

#dimer
def antibody_antigen(meta):
    res = QueryAntibody(params['verbose']).combo2_abag()
    for v in res.values():
        rec = (v['pdb_id'], v['combo_id'], v['chain_combo'], 
            v['chain_label'], 2, 'AbAg')
        meta['antibody-antigen'] += 1
        yield rec

#trimer
def heavy_light_antigen(meta):
    res = QueryAntibody(params['verbose']).combo3_heavy_light_antigen()
    for v in res.values():
        rec = (v['pdb_id'], v['combo_id'], v['chain_combo'], 
            v['chain_label'], 3, 'HLAg')
        meta['heavy-light-antigen'] += 1
        yield rec


if __name__ == "__main__":
    params = {
        'verbose': True,
        'table_name': 'pdb_complex',
    }
    UpdateComplex(params['table_name'], params['verbose']).empty_tables()

    pool = [
        # antibody complex
        heav_light, antibody_antigen, heavy_light_antigen,
    ]
    for func in pool:
        record_iter = func(meta)
        bc = BuildComplex(record_iter, params['verbose'])
        bc.insert_batch_records(params['table_name'])

    footer(params, meta)
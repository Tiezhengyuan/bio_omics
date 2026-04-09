#!/usr/bin/python
'''
requirements:
    - meta_pdb
    - table view_antibody
example:
    - parallel -j16 python app.py pdb_complex_antibody ::: {0..237}
    - python app.py pdb_complex_antibody 0
function:
    - identify complex
    - insert data into table pdb_complex
'''
from src.helper import *
from bioomics import QueryAntibody

def sort_ids(item:dict):
    df = pd.DataFrame({
        'no': item['chain_combo'],
        'label': item['chain_label'],
        'ids': item['chain_ids'],
    }).sort_values('no')
    pdb_id = item['pdb_id']
    model_no = int(item['model_no'])
    chain_combo = df['no'].str.cat(sep='-')
    chain_label = df['label'].str.cat(sep='-')
    combo_id = f"{pdb_id}_{item['model_no']}_{chain_combo}"
    combo = len(df)
    chain_ids = df['ids'].str.cat(sep=',')
    return pdb_id, model_no, combo_id, chain_combo, chain_label, combo, chain_ids


def heav_light(params, meta):
    category = 'AbAb'
    combo_ids = QueryComplex(params['verbose']).combo_ids(params['pdb_group'], category)

    qa = QueryAntibody(params['verbose'])
    res = qa.dimer_heavy_light(params['pdb_group'])
    for item in res:
        rec = sort_ids(item) + (category,)
        if rec[2] not in combo_ids:
            meta['heavy-light'] += 1
            yield rec

def chain_target(params, meta):
    category = 'AbAg'
    combo_ids = QueryComplex(params['verbose']).combo_ids(params['pdb_group'], category)

    qa = QueryAntibody(params['verbose'])
    res = qa.dimer_antibody_target(params['pdb_group'])
    for item in res:
        rec = sort_ids(item)+ (category,)
        if rec[2] not in combo_ids:
            meta['chain-target'] += 1
            yield rec

def heavy_light_target(params, meta):
    category = 'Ab2Ag'
    combo_ids = QueryComplex(params['verbose']).combo_ids(params['pdb_group'], category)

    qa = QueryAntibody(params['verbose'])
    res = qa.trimer_heavy_light_target(params['pdb_group'])
    for item in res:
        rec = sort_ids(item) + (category,)
        if rec[2] not in combo_ids:
            meta['heavy-light-target'] += 1
            yield rec

def heavy2_light2(params, meta):
    category = 'Ab2Ab2'
    combo_ids = QueryComplex(params['verbose']).combo_ids(params['pdb_group'], category)

    qa = QueryAntibody(params['verbose'])
    res = qa.tetramer_antibody(params['pdb_group'])
    for item in res:
        rec = sort_ids(item) + (category,)
        if rec[2] not in combo_ids:
            meta['antibody'] += 1
            yield rec

def antibody_target(params, meta):
    category = 'Ab2Ab2Ag'
    combo_ids = QueryComplex(params['verbose']).combo_ids(params['pdb_group'], category)

    qa = QueryAntibody(params['verbose'])
    res = qa.pentamer_antibody_target(params['pdb_group'])
    for item in res:
        rec = sort_ids(item) + (category,)
        if rec[2] not in combo_ids:
            meta['antibody-target'] += 1
            yield rec

if __name__ == "__main__":
    params.update({
        'verbose': True,
        'table_name': 'pdb_complex',
    })

    pool = [
        heav_light, 
        chain_target,
        heavy_light_target,
        heavy2_light2,
        antibody_target,
    ]
    for func in pool:
        record_iter = func(params, meta)
        bc = BuildComplex(record_iter, params['verbose'])
        bc.insert_batch_records(params['table_name'])

    footer(params, meta)
#!/usr/bin/python
'''
requirements:
    - view_antibody
example:
    - parallel -j8 bomics label_antibody ::: {0..241}
    - bomics label_antibody 0
function:
    - identify protein and insert data into table label_chain
    - identify complex and insert data into table label_complex
'''
from src.helper import *
from bioomics import QueryAntibody

def ab_chain(params, meta):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        SELECT pdb_id, chain_id, first_chain_id, chain_type
        FROM view_antibody
        WHERE pdb_group = {pdb_group}
            AND chain_id NOT IN (
                SELECT chain_id from {table_name}
            )
    """
    rows = QueryComplex(params['verbose']).list_data(query)
    meta['ab_chain'] += len(rows)
    for row in rows:
        rec = (row['pdb_id'], row['chain_id'], row['first_chain_id'],
            row['chain_type'], 'Ab',)
        yield rec

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
    qc = QueryComplex(params['verbose'])
    combo_ids = qc.combo_ids(params['pdb_group'], category)

    qa = QueryAntibody(params['verbose'])
    res = qa.dimer_heavy_light(params['pdb_group'])
    for item in res:
        rec = sort_ids(item) + (category,)
        if rec[2] not in combo_ids:
            meta['heavy-light'] += 1
            yield rec

def chain_target(params, meta):
    category = 'AbAg'
    qc = QueryComplex(params['verbose'])
    combo_ids = qc.combo_ids(params['pdb_group'], category)

    qa = QueryAntibody(params['verbose'])
    res = qa.dimer_antibody_target(params['pdb_group'])
    for item in res:
        rec = sort_ids(item)+ (category,)
        if rec[2] not in combo_ids:
            meta['chain-target'] += 1
            yield rec

def heavy_light_target(params, meta):
    category = 'Ab2Ag'
    qc = QueryComplex(params['verbose'])
    combo_ids = qc.combo_ids(params['pdb_group'], category)

    qa = QueryAntibody(params['verbose'])
    res = qa.trimer_heavy_light_target(params['pdb_group'])
    for item in res:
        rec = sort_ids(item) + (category,)
        if rec[2] not in combo_ids:
            meta['heavy-light-target'] += 1
            yield rec

def heavy2_light2(params, meta):
    category = 'Ab2Ab2'
    qc = QueryComplex(params['verbose'])
    combo_ids = qc.combo_ids(params['pdb_group'], category)

    qa = QueryAntibody(params['verbose'])
    res = qa.tetramer_antibody(params['pdb_group'])
    for item in res:
        rec = sort_ids(item) + (category,)
        if rec[2] not in combo_ids:
            meta['antibody'] += 1
            yield rec

def antibody_target(params, meta):
    category = 'Ab2Ab2Ag'
    qc = QueryComplex(params['verbose'])
    combo_ids = qc.combo_ids(params['pdb_group'], category)

    qa = QueryAntibody(params['verbose'])
    res = qa.pentamer_antibody_target(params['pdb_group'])
    for item in res:
        rec = sort_ids(item) + (category,)
        if rec[2] not in combo_ids:
            meta['antibody-target'] += 1
            yield rec

if __name__ == "__main__":
    params.update({
        'chunk_size': 5,
    })

    # label chain
    pool = [
        ab_chain,
    ]
    for func in pool:
        params['table_name'] = 'label_chain'
        record_iter = func(params, meta)
        bc = BuildComplex(params['verbose'], params['chunk_size'])
        bc.insert_batch_records(record_iter, params['table_name'])

    # label complex
    pool = [
        heav_light, 
        chain_target,
        heavy_light_target,
        heavy2_light2,
        antibody_target,
    ]
    for func in pool:
        params['table_name'] = 'label_complex'
        record_iter = func(params, meta)
        bc = BuildComplex(params['verbose'], params['chunk_size'])
        bc.insert_batch_records(record_iter, params['table_name'])

    footer(params, meta)
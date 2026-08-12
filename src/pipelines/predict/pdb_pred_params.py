#!/usr/bin/python
'''
requirements:
    - table label_chain
example:
    - parallel -j8 bomics pdb_pred_params ::: {0..241}
    - bomics pdb_pred_params 0
function:
    - build params for pdb structure prediction
'''
from src.helper import *
from .func_predict import *

def build(predictor_name, rows, params, meta):
    for row in rows:
        pdb_id = row['pdb_id']
        query_faa = params['input_dir'] / row['relative_faa']
        if query_faa.is_file():
            meta[f"{params['category']}_records"] += 1
            outdir = params['outdir'] / pdb_id[:2] /pdb_id / row['input_id']
            rec = (
                pdb_id,
                row['input_id'],
                params['input_table_name'],
                params['input_col_name'],
                params['category'],
                predictor_name,
                str(query_faa),
                str(outdir),
            )
            yield rec
        else:
            meta[f"{params['category']}_skip"] += 1

if __name__ == "__main__":
    params.update({
        'table_name': 'pdb_pred_params',
    })
    add_params(params)

    # AF2
    qc = QueryComplex(params['verbose'])
    pool = (
        retrieve_chain(params, 'chain_af'),
        retrieve_combo(params, 'combo2_af', 'combo2'),
        retrieve_combo(params, 'combo3_af', 'combo3'),
        retrieve_combo(params, 'combo4_af', 'combo4'),
        retrieve_combo(params, 'combo5_af', 'combo5'),
    )
    for curr_params, query in pool:
        rows = qc.list_data(query)
        print(f"{len(rows)} PDB files ready...")
        meta[f"{curr_params['outdir_name']}_rows"] += len(rows)        

        record_iter = build('AF2', rows, curr_params, meta)
        bc = BuildComplex(params['verbose'])
        bc.insert_batch_records(record_iter, params['table_name'])
    
    footer(meta, params)
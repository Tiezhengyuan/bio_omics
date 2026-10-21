#!/usr/bin/python
'''
requirements:
    - table af_msa
example:
    - parallel -j8 bomics af_pred ::: {0..241}
    - bomics af_pred 0
function:
    - prefix monomer
    - insert data into table af_pred
'''
from .func_predict import *


if __name__ == "__main__":
    params.update({
        'tool_name': 'colabfold_batch',
        'num_recycle': '3',
        'relax': True,
        'max_iterations': '2000',
        'table_name': 'chain_af',
    })
    add_params(params)
    validate_predict_data(params, meta)

    rows = retrieve_chain(params, meta)
    for row in rows:
        # params for this run
        curr_params = dict(params)
        chain_params(row, curr_params)
        
        # run predictor
        run_colabfold_batch(curr_params, meta)

    footer(meta, params)
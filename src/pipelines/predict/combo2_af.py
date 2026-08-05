#!/usr/bin/python
'''
requirements:
    - combo2
example:
    - parallel -j8 bomics combo2_af ::: {0..241}
    - bomics combo2_af 0
function:
    - predict dimer complex
    - insert data into table combo2_af
'''

from src.helper import *
from .func_predict import *

if __name__ == "__main__":
    params.update({
        'tmp_dir': os.getenv('tmp_dir'),
        'tool_name': 'colabfold_batch',
        'combo': 2,
        'table_name': 'combo2_af',
        'num_recycle': '3',
        'relax': True,
    })
    add_params(params)
    # remove rows where no exports
    validate_predict_data(params, meta)

    rows = retrieve_combo(params, meta)
    for row in rows:
        # build params only for this run
        curr_params = dict(params)
        combo_params(row, curr_params)
        # predict
        run_colabfold_batch(curr_params, meta)

    footer(meta, params)
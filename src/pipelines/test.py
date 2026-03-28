#!/usr/bin/python
'''
requirements:
    - meta_combo2
example:
    - parallel --delay 5 -j 8 python app.py test ::: {0..117}
    - python app.py test
function:
    - test concurrency
'''

import time
from src.helper import *
from bioomics import Predictor


if __name__ == "__main__":
    args = {
        'script': sys.argv[0],
        'pdb_group': None,
        'overwrite': False,
        'verbose': True,
        'test': False,
        'data_dir': os.getenv('data_dir'),
        'predict_dir': os.getenv('predict_dir'),
        'table_name': 'combo2_alphafold',
    }
    add_meta(args)


    #update record
    time.sleep(30)
    UpdateComplex('meta_execute', args['verbose']).update_execute_time(args)

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

if __name__ == "__main__":
    params = {
        'script': sys.argv[0],
        'pdb_group': None,
        'overwrite': False,
        'verbose': True,
        'test': False,
        'data_dir': os.getenv('data_dir'),
        'predict_dir': os.getenv('predict_dir'),
        'table_name': 'meta',
    }
    add_params(params)


    #update record
    time.sleep(30)
    footer(params, meta)

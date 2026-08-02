#!/usr/bin/python
'''
requirements:
    - table combo4
example:
    - parallel -j8 bomics combo4_relaxed ::: {0..241}
    - bomics combo4_relaxed 0
function:
    - fix pdb files
    - relax pdb
    - insert data into table combo4_relaxed
'''
from src.helper import *
from .func_relaxed import *

if __name__ == "__main__":
    params.update({
        'dry_run': False,
        'overwrite': False,
        'tmp_dir': os.getenv('tmp_dir'),
        'combo': 4,
        'table_name': 'combo4_relaxed',
        'force_fields': ['amber14-all.xml', 'amber14/tip3pfb.xml'],
        'max_iterations': 2_000,
        'steps': 10_000,
    })
    
    rows = retrieve_pdb(params, meta)
    if params['dry_run']:
        record_iter = dry_run(rows, params, meta)
        bc = BuildComplex(params['verbose'], 2)
        bc.insert_batch_records(record_iter, params['table_name'])
    else:
        for row in rows:
            rec = fix_and_relax(row, params, meta)
            bc = BuildComplex(params['verbose'])
            bc.insert_batch_records([rec,], params['table_name'])

    footer(meta)
'''
requirements:
    - The table label_complex is ready.
    - meta_chain_pdb is ready
example:
    - parallel -j8 bomics combo4 ::: {0..241}
    - bomics combo4 0
function:
    - split PDB of complex to paired chains PDB.
    - export complex csv and fasta for structure predication
note:
    - only consider model.id == 0
'''
from src.helper import *
from .func_combo import *

##################################
if __name__ == "__main__":
    params.update({
        'overwrite': False,
        'combo': 4,
        'table_name': 'combo4',
    })

    row_iter = retrieve_data(params, meta)
    record_iter = build_data(row_iter, params, meta)
    
    bc = BuildComplex(params['verbose'])
    bc.insert_batch_records(record_iter, params['table_name'])

    footer(meta)
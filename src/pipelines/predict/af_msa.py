#!/usr/bin/python
'''
requirements:
    - table label_chain
example:
    - parallel -j8 bomics af_msa ::: {0..241}
    - bomics af_msa 0
function:
    - MSA and build *.a3m for structure prediction
'''
import asyncio
from bioomics import async_concurr
from .func_predict import *

def retrieve(params, meta):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        SELECT p.pdb_id, p.input_id, p.query_faa, p.outdir
        FROM pdb_pred_params p
        LEFT JOIN meta_pdb   m ON p.pdb_id=m.pdb_id
        WHERE m.pdb_group={pdb_group}
            AND p.predictor = 'AF2'
            AND NOT EXISTS (
                SELECT 1 FROM {table_name} b
                WHERE p.input_id = b.input_id
            )
        ORDER BY p.category
    """
    rows = QueryComplex(params['verbose']).list_data(query)
    print(f"{len(rows)} query files ready...")
    meta["rows"] += len(rows)
    wrapped = []
    for row in rows:
        row['a3m_name'] = f"{row['input_id']}.a3m"
        row['a3m_file'] = Path(row['outdir']) / row['a3m_name']
        wrapped.append(row)
    return wrapped

# NO asynchronous funciton is needed
def post_func(input):
    # input format is consistent with returns of 
    # run_colabfold_search and run_colabfold_msa
    row, params, meta = input
    if row.get('a3m_file') and row['a3m_file'].is_file():
        # insert
        rec = (row['input_id'], row['a3m_name'], params['host'],)
        if params['verbose']:
            print('### Try to insert ', rec)
        bc = BuildComplex(params['verbose'])
        bc.insert_batch_records([tuple(rec),], params['table_name'])
        meta['records'] += 1
    else:
        meta['no_a3m'] += 1


if __name__ == "__main__":
    params.update({
        'tmp_dir': os.getenv('tmp_dir'),
        'colabfold_db': os.getenv('colabfold_db'),
        'overwrite':False,
        'params_table_name': 'pdb_pred_params',
        'table_name': 'af_msa',
    })
    add_params(params)
    validate_file_rows('a3m_name', params, meta)
    
    # asynchronous functions
    funcs = [
        run_colabfold_search,
        run_colabfold_msa,
    ]
    rows = retrieve(params, meta)
    input = [(row, params, meta) for row in rows]
    asyncio.run(
        async_concurr.main(funcs, input, post_func)
    )
    
    footer(meta, params)
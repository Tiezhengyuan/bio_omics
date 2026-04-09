'''
requirements:
    - bfactor.pkl
    - table meta_bfactor
example: 
    - parallel -j4 python app.py pdb_bfactor ::: {0..237}
    - python app.py pdb_bfactor 0
functions:
    - put average bfactor to table pdb_bfactor
'''

from src.helper import *

@header
def retrieve_pdb(params, meta):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        SELECT B.pdb_id, B.relative_pkl
        FROM meta_bfactor B, meta_pdb P
        WHERE B.pdb_id=P.pdb_id
            AND P.pdb_group = {pdb_group}
            AND B.relative_pkl IS NOT NULL
            AND B.pdb_id NOT IN (
                SELECT DISTINCT pdb_id FROM {table_name}
            );
    """
    # each record must have pdb_id and pdb_file
    rows = QueryComplex(params['verbose']).list_data(query)
    print(f"{len(rows)} files ready...")
    meta['pdb_group'] = int(pdb_group)
    meta['rows'] += len(rows)
    for row in rows:
        row['pickle_file'] = os.path.join(params['data_dir'], row['relative_pkl'])
        if os.path.isfile(row['pickle_file']):
            yield row

if __name__ == "__main__":
    params.update({
        'table_name': 'pdb_bfactor',
    })
    params_iter = retrieve_pdb(params, meta)

    # put bfactor into table
    print(f"Try put meta data to {params['table_name']}...")
    bc = BuildComplex(params_iter, params['verbose'])
    bc.insert_pdb_bfactor()

    footer(meta)
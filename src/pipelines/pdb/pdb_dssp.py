'''
requirements:
    - PDB raw data
    - table pdb, meta_chain_pdb
example: 
    - parallel -j4 python app.py pdb_dssp ::: {0..237}
    - python app.py pdb_dssp 0
functions:
    - Analyze secondary structure by chains
    - put dssp sequence to table pdb_dssp
'''

from src.helper import *

@header
def search_pdb(params, meta):
    pdb_group = params['pdb_group']
    # 1. ==pdb_group, 2. pdb_id not in pdb_secondary 3. chain is protein
    table_name = params['table_name']
    query = f"""
        SELECT I.pdb_id, I.model_no, I.chain_no, P.chain_pdb
        FROM pdb_chainid I
        LEFT JOIN pdb_chain_seq S  ON I.chain_id = S.chain_id
        LEFT JOIN meta_chain_pdb P ON I.chain_id = P.chain_id
        LEFT JOIN meta_pdb M       ON I.pdb_id = M.pdb_id
        WHERE M.pdb_group = {pdb_group}
            AND P.chain_pdb IS NOT NULL
            AND S.chain_seq IS NOT NULL
            AND I.model_no = 0
            AND I.pdb_id NOT IN (
                SELECT pdb_id FROM {table_name}
            )
    ;"""
    rows = QueryComplex(params['verbose']).list_data(query)
    meta['rows'] = len(rows)
    print(f"Group {pdb_group}: Collect chains from {len(rows)} PDB ...")
    for row in rows:
        pdb_file = os.path.join(params['data_dir'], row['chain_pdb'])
        if os.path.isfile(pdb_file):
            row['pdb_file'] = pdb_file
            yield row

if __name__ == "__main__":
    params.update({
        'chunk_size': 10,
        'table_name': 'pdb_dssp',
    })
    row_iter = search_pdb(params, meta)

    # put bfactor into table pdb_bfactor
    print(f"Try to analyze secondary structure and insert results into {params['table_name']}...")
    bc = BuildComplex(row_iter, params['verbose'], params['chunk_size'])
    bc.insert_pdb_dssp()

    footer(meta)
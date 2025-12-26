'''
requirements:
    - PDB raw data
    - table pdb, meta_chain_pdb
example: 
    - parallel -j4 python src/pipelines/pdb_dssp.py ::: {0..116}
    - python src/pipelines/pdb_dssp.py 0
functions:
    - Analyze secondary structure by chains
    - put dssp sequence to table pdb_dssp
'''

from helper import *

def search_pdb(pdb_group, args, meta):
    # 1. ==pdb_group, 2. pdb_id not in pdb_secondary 3. chain is protein
    table_name = args['table_name']
    query = f"""
        SELECT I.pdb_id, I.model_no, I.chain_no, P.relative_path
        FROM pdb_chainid I
        LEFT JOIN pdb_chain_seq S  ON I.chain_id = S.chain_id
        LEFT JOIN meta_chain_pdb P ON I.chain_id = P.chain_id
        LEFT JOIN meta_pdb M       ON I.pdb_id = M.pdb_id
        WHERE M.pdb_group = {pdb_group}
            AND P.relative_path IS NOT NULL
            AND S.chain_seq IS NOT NULL
            AND I.model_no = 0
            AND I.pdb_id NOT IN (
                SELECT pdb_id FROM {table_name}
            )
    ;"""
    rows = QueryComplex(args['verbose']).list_data(query)
    meta['pdb_group'] = pdb_group
    meta['rows'] = len(rows)
    print(f"Group {pdb_group}: Collect chains from {len(rows)} PDB ...")
    for row in rows:
        pdb_file = os.path.join(args['data_dir'], row['relative_path'])
        if os.path.isfile(pdb_file):
            row['pdb_file'] = pdb_file
            yield row

if __name__ == "__main__":
    args = {
        'verbose': True,
        'chunk_size': 10,
        'data_dir': os.getenv('data_dir'),
        'table_name': 'pdb_dssp',
    }
    pdb_group = sys.argv[1]
    row_iter = search_pdb(pdb_group, args, meta)

    # put bfactor into table pdb_bfactor
    print(f"Try to analyze secondary structure and insert results into {args['table_name']}...")
    bc = BuildComplex(row_iter, args['verbose'], args['chunk_size'])
    bc.insert_pdb_dssp()

    print(meta)
'''
requirements:
    - bfactor.pkl
    - table meta_bfactor
example: 
    - parallel -j4 python src/pipelines/pdb_bfactor.py ::: {0..116}
    - python src/pipelines/pdb_bfactor.py 0
functions:
    - put average bfactor to table pdb_bfactor
'''

from helper import *

def retrieve_pdb(pdb_group, args, meta):
    table_name = args['table_name']
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
    rows = QueryComplex(args['verbose']).list_data(query)
    print(f"{len(rows)} files ready...")
    meta['pdb_group'] = int(pdb_group)
    meta['rows'] += len(rows)
    for row in rows:
        row['pickle_file'] = os.path.join(args['data_dir'], row['relative_pkl'])
        if os.path.isfile(row['pickle_file']):
            yield row


if __name__ == "__main__":
    args = {
        'verbose': True,
        'data_dir': os.getenv('data_dir'),
        'table_name': 'pdb_bfactor',
    }
    pdb_group = sys.argv[1]
    args_iter = retrieve_pdb(pdb_group, args, meta)

    # put bfactor into table
    print(f"Try put meta data to {args['table_name']}...")
    bc = BuildComplex(args_iter, args['verbose'])
    bc.insert_pdb_bfactor()

    print(meta)


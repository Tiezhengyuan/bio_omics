'''
requirements:
    - blastp alignment *.blastp5 files
    - table meta_blastp5
example: 
    - parallel -j4 python app.py align_uniprot ::: {0..237}
    - python app.py align_uniprot 0
functions:
    - pust alignment results of blastp into align_uniprot
'''

from src.helper import *

@header
def retrieve_pdb(args, meta):
    pdb_group = args['pdb_group']
    table_name = args['table_name']
    query = f"""
        SELECT B.pdb_id, B.relative_aln
        FROM meta_blastp5 B
        LEFT JOIN meta_pdb P ON B.pdb_id=P.pdb_id
        WHERE P.pdb_group = {pdb_group}
            AND B.pdb_id NOT IN (
                SELECT DISTINCT pdb_id FROM {table_name}
           )
    ;"""
    rows = QueryComplex(args['verbose']).list_data(query)
    meta['align_uniprot'] = len(rows)
    for row in rows:
        row['aln_file'] = os.path.join(args['data_dir'], row['relative_aln'])
        yield row

if __name__ == "__main__":
    args = {
        'pdb_group': sys.argv[1],
        'verbose': True if os.getenv('verbose') == 'true' else False,
        'data_dir': os.getenv('data_dir'),
        'table_name': 'align_uniprot',
    }
    record_iter = retrieve_pdb(args, meta)

    print(f"Try to insert alignment into table {args['table_name']}")
    bc = BuildComplex(record_iter, args['verbose'])
    bc.insert_align_uniprot()

    footer(meta)

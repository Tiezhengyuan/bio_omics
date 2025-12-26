'''
requirements:
    - blastp alignment *.blastp5 files
    - table meta_blastp5
example: 
    - parallel -j4 python src/pipelines/align_uniprot.py ::: {0..116}
    - python src/pipelines/align_uniprot.py 0
functions:
    - pust alignment results of blastp into align_uniprot
'''

from helper import *

def retrieve_pdb(pdb_group, args, meta):
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
    meta['pdb_group'] = int(pdb_group)
    meta['align_uniprot'] = len(rows)
    for row in rows:
        row['aln_file'] = os.path.join(args['data_dir'], row['relative_aln'])
        yield row

if __name__ == "__main__":
    args = {
        'verbose': True,
        'data_dir': os.getenv('data_dir'),
        'table_name': 'align_uniprot',
    }
    pdb_group = sys.argv[1]
    record_iter = retrieve_pdb(pdb_group, args, meta)

    print(f"Try to insert alignment into table {args['table_name']}")
    bc = BuildComplex(record_iter, args['verbose'])
    bc.insert_align_uniprot()

    print(meta)

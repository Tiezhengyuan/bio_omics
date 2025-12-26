#!/usr/bin/python
'''
requirements:
    - table meta_faa
    - igblastp ready
example:
    - parallel -j4 python src/pipelines/align_vregion.py ::: {0..116}
    - python src/pipelines/align_vregion.py 0
function:
    - build table pdb and insert data into table align_vregion
'''

from helper import *

def retrieve_data(pdb_group, args, meta):
    table_name = args['table_name']
    query = f"""
        WITH tmp AS (
            SELECT F.pdb_id, F.relative_faa, I.relative_aln 
            FROM pdb_faa F, meta_igblastp I
            WHERE F.pdb_id = I.pdb_id
                AND F.model_no = I.model_no
                AND F.relative_faa IS NOT NULL
                AND I.relative_aln IS NOT NULL
        )
        SELECT T.pdb_id, T.relative_faa, T.relative_aln
        FROM tmp T
        LEFT JOIN meta_pdb P ON T.pdb_id = P.pdb_id
        WHERE P.pdb_group = {pdb_group}
            AND T.pdb_id NOT IN(
                SELECT DISTINCT pdb_id FROM {table_name}
            )
    """
    rows = QueryComplex(args['verbose']).list_data(query)
    meta['pdb_group'] = int(pdb_group)
    meta['rows'] = len(rows)
    for row in rows:
        row['faa_file'] = os.path.join(args['data_dir'], row['relative_faa'])
        row['aln_file'] = os.path.join(args['data_dir'], row['relative_aln'])
        yield row

if __name__ == "__main__":
    args = {
        'verbose': True,
        'data_dir': os.getenv('data_dir'),
        'table_name': 'align_vregion',
    }
    pdb_group = sys.argv[1]
    row_iter = retrieve_data(pdb_group, args, meta)

    print(f"Try to insert data into {args['table_name']}")
    bc = BuildComplex(row_iter, args['verbose'])
    bc.insert_align_vregion()

    print(meta)


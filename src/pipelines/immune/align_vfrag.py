#!/usr/bin/python
'''
requirements:
    - table meta_faa
    - igblastp ready
example:
    - parallel -j4 python app.py align_vfrag ::: {0..237}
    - python app.py align_vfrag 0
function:
    - build table pdb and insert FRs and CDRs of V-region into table align_vfrag
'''

from src.helper import *

@header
def retrieve(params, meta):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        WITH tmp AS (
            SELECT F.pdb_id, F.relative_faa, I.v_aln
            FROM pdb_faa F, meta_vdj I
            WHERE F.pdb_id = I.pdb_id
                AND F.model_no = I.model_no
                AND F.relative_faa IS NOT NULL
        )
        SELECT T.pdb_id, T.relative_faa, T.v_aln
        FROM tmp T
        LEFT JOIN meta_pdb P ON T.pdb_id = P.pdb_id
        WHERE P.pdb_group = {pdb_group}
            AND T.pdb_id NOT IN(
                SELECT DISTINCT pdb_id FROM {table_name}
            )
    """
    rows = QueryComplex(params['verbose']).list_data(query)
    meta['rows'] = len(rows)
    for row in rows:
        row['faa_file'] = os.path.join(params['data_dir'], row['relative_faa'])
        row['aln_file'] = os.path.join(params['data_dir'], row['v_aln'])
        if os.path.isfile(row['faa_file']) and os.path.isfile(row['aln_file']):
            yield row


if __name__ == "__main__":
    params.update({
        'table_name': 'align_vfrag',
    })
    print(f"Try to insert data into {params['table_name']}...")
    row_iter = retrieve(params, meta)
    bc = BuildComplex(row_iter, params['verbose'])
    bc.insert_align_vfrag()


    footer(meta)


#!/usr/bin/python
'''
requirements:
    - table meta_faa
    - igblastp ready
example:
    - parallel -j4 python app.py align_dregion ::: {0..117}
    - python app.py align_dregion 0
function:
    - build table pdb and insert D-region into table align_dregion
'''

from src.helper import *
from bioomics import ParseBlastXml

@header
def retrieve(params, meta):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        WITH tmp AS (
            SELECT F.pdb_id, F.relative_faa, I.d_aln
            FROM pdb_faa F, meta_vdj I
            WHERE F.pdb_id = I.pdb_id
                AND F.model_no = I.model_no
                AND F.relative_faa IS NOT NULL
        )
        SELECT T.pdb_id, T.relative_faa, T.d_aln
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
        row['aln_file'] = os.path.join(params['data_dir'], row['d_aln'])
        if os.path.isfile(row['faa_file']) and os.path.isfile(row['aln_file']):
            yield row

def run(row_iter, meta):
    for row in row_iter:
        hits = ParseBlastXml(row['aln_file']).hits()
        for hit in hits:
            hit_id, organism = None, None
            if 'hit_id' in hit:
                items = hit['hit_id'].split('|')
                hit_id = items[0]
                if len(items) > 1:
                    organism = items[1]
            rec = (row['pdb_id'], hit['query'], hit_id, organism,
                hit.get('evalue'),
                hit.get('bit_score'),
                hit.get('identities'),
                hit.get('query_len'),
                hit.get('align_len'),
                hit.get('query_start'),
                hit.get('query_end'),
                hit.get('query_seq'),
                hit.get('match_seq'),
                hit.get('sbjct_seq'),
                hit.get('sbjct_start'),
                hit.get('sbjct_end'),
            )
            meta['hits'] += 1
            yield rec

if __name__ == "__main__":
    params = {
        'pdb_group': sys.argv[1],
        'verbose': True,
        'data_dir': os.getenv('data_dir'),
        'table_name': 'align_dregion',
    }
    print(f"Try to insert data into {params['table_name']}...")
    row_iter = retrieve(params, meta)
    record_iter = run(row_iter, meta)
    bc = BuildComplex(record_iter, params['verbose'])
    bc.insert_batch_records(params['table_name'])

    print(meta)


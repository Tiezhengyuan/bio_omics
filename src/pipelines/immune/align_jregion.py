#!/usr/bin/python
'''
requirements:
    - table meta_faa
    - igblastp ready
example:
    - parallel -j4 python app.py align_jregion ::: {0..117}
    - python app.py align_jregion 0
function:
    - build table pdb and insert J-region into table align_jregion
'''

from src.helper import *
from bioomics import ParseBlastXml

@header
def retrieve(params, meta):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        WITH tmp AS (
            SELECT F.pdb_id, F.relative_faa, I.j_aln
            FROM pdb_faa F, meta_vdj I
            WHERE F.pdb_id = I.pdb_id
                AND F.model_no = I.model_no
                AND F.relative_faa IS NOT NULL
        )
        SELECT T.pdb_id, T.relative_faa, T.j_aln
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
        row['aln_file'] = os.path.join(params['data_dir'], row['j_aln'])
        if os.path.isfile(row['faa_file']) and os.path.isfile(row['aln_file']):
            yield row

def run(row_iter, meta):
    for row in row_iter:
        hits = ParseBlastXml(row['aln_file']).hits()
        for hit in hits:
            meta['hits'] += 1
            item = {
                'pdb_id': row['pdb_id'],
                'chain_id': hit['query'],
                'evalue': hit.get('evalue'),
                'bit_score': hit.get('bit_score'),
                'identities': hit.get('identities'),
                'query_len': hit.get('query_len'),
                'align_len': hit.get('align_len'),
                'query_start': hit.get('query_start'),
                'query_end': hit.get('query_end'),
                'query_seq': hit.get('query_seq'),
                'match_seq': hit.get('match_seq'),
                'sbjct_seq': hit.get('sbjct_seq'),
                'sbjct_start': hit.get('sbjct_start'),
                'sbjct_end': hit.get('sbjct_end'),
            }
            item.update(
                ParseBlastXml.parse_imgt_allele(hit['hit_id'])
            )
            rec = [item[col] for col in params['table_cols']]
            yield tuple(rec)


if __name__ == "__main__":
    params.update({
        'table_name': 'align_jregion',
        'table_cols': ['pdb_id', 'chain_id', 'allele_name', 'specie',
            'gene_name', 'gene_family', 'chain_type', 'isotype',
            'evalue', 'bit_score', 'identities', 'query_len', 'align_len',
            'query_start', 'query_end', 'query_seq', 'match_seq',
            'sbjct_seq', 'sbjct_start', 'sbjct_end',],
    })
    print(f"Try to insert data into {params['table_name']}...")
    row_iter = retrieve(params, meta)
    record_iter = run(row_iter, meta)
    bc = BuildComplex(record_iter, params['verbose'])
    bc.insert_batch_records(params['table_name'], params['table_cols'])

    footer(meta)


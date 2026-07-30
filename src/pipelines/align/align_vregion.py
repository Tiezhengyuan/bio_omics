#!/usr/bin/python
'''
requirements:
    - table meta_faa
    - igblastp ready
example:
    - parallel -j8 bomics align_vregion ::: {0..241}
    - bomics align_vregion 0
function:
    - build table pdb and insert data into table align_vregion
'''

from src.helper import *
from bioomics import ParseIgblast, TABLE_COLUMNS

@header
def retrieve_data(params, meta):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        WITH tmp AS (
            SELECT F.pdb_id, F.relative_faa, I.v_aln 
            FROM  pdb_faa F, align_vdj I
            WHERE F.pdb_id = I.pdb_id
                AND F.model_no = I.model_no
                AND F.relative_faa IS NOT NULL
                AND I.v_aln IS NOT NULL
        )
        SELECT T.pdb_id, T.relative_faa, T.v_aln
        FROM tmp           T
        LEFT JOIN meta_pdb P ON T.pdb_id = P.pdb_id
        WHERE P.pdb_group = {pdb_group}
            AND T.pdb_id NOT IN(
                SELECT DISTINCT pdb_id FROM {table_name}
            )
    ;"""
    rows = QueryComplex(params['verbose']).list_data(query)
    meta['rows'] = len(rows)
    for row in rows:
        row['faa_file'] = params['output_pdb_dir'] / row['relative_faa']
        row['aln_file'] = params['output_protein_dir'] / row['v_aln']
        if row['faa_file'].is_file() and row['aln_file'].is_file():
            yield row


def align_vregion(row_iter, params, meta):
    cols = TABLE_COLUMNS[params['table_name']]
    for row in row_iter:
        pi = ParseIgblast(row['aln_file'], row['faa_file'], params['verbose'])
        for hit in pi.retrieve_hits():
            item = {
                'pdb_id': row['pdb_id'],
                'chain_id': hit['query_id'],
            }
            if 'allele_name' in hit:
                allele_name, specie = hit['allele_name'].split('|')
                region_names = [i.replace('-IMGT', '') for i in hit.get('region_names', [])]
                item.update({
                    'allele_name': allele_name,
                    'specie': specie,
                    'gene_name': hit['gene_name'],
                    'gene_family': hit['gene_family'],
                    'chain_type': hit['chain_type'].upper(),
                    'isotype': hit['isotype'].upper(), 
                    'evalue': float(hit['e-value']),
                    'bit_score': float(hit['bits-score']),
                    'region_names': ','.join(region_names) if region_names else None,
                })
                meta['records'] += 1
            else:
                meta['null'] += 1
            rec = [item.get(col) for col in cols]
            if params['verbose']:
                print(rec)
            yield tuple(rec)

if __name__ == "__main__":
    params.update({
        'chunk_size': 50,
        'table_name': 'align_vregion',
    })

    row_iter = retrieve_data(params, meta)
    record_iter = align_vregion(row_iter, params, meta)

    print(f"Try to insert data into {params['table_name']}")
    bc = BuildComplex(params['verbose'], params['chunk_size'])
    bc.insert_batch_records(record_iter, params['table_name'])

    footer(meta)
#!/usr/bin/python
'''
requirements:
    - table meta_faa
    - igblastp ready
example:
    - parallel -j8 python app.py align_vfrag ::: {0..240}
    - python app.py align_vfrag 0
function:
    - build table pdb and insert FRs and CDRs of V-region into table align_vfrag
'''

from src.helper import *
from bioomics import ParseIgblast, TABLE_COLUMNS

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
        row['faa_file'] = os.path.join(params['output_pdb_dir'], row['relative_faa'])
        row['aln_file'] = os.path.join(params['output_pdb_dir'], row['v_aln'])
        if os.path.isfile(row['faa_file']) and os.path.isfile(row['aln_file']):
            yield row

def align_vfrag(row_iter, params, meta):
    cols = TABLE_COLUMNS[params['table_name']]
    for row in row_iter:
        # row: pdb_id, aln_file, faa_file
        c = ParseIgblast(row['aln_file'], row['faa_file'], params['verbose'])
        for region in c.retrieve_regions():
            item = {
                'pdb_id': row['pdb_id'],
                'chain_id': region['query_id'],
            }
            if 'percent identity' in region:
                region_name, seq = region['region'], region['seq']
                item['identity'] = region['percent identity']
                item['seq'] = seq
                if  region_name == 'Total':
                    item.update({
                        'fragment': 'Total',
                        'seq_from': 1,
                        'seq_to': len(seq),
                    })
                else:
                    item.update({
                        'fragment': region_name.replace('-IMGT', ''),
                        'seq_from': region['from'],
                        'seq_to': region['to'],
                    })
                meta['fragment'] += 1
            # build record
            rec = [item.get(col) for col in cols]
            if params['verbose']:
                print(rec)
            yield tuple(rec)


if __name__ == "__main__":
    params.update({
        'table_name': 'align_vfrag',
    })
    row_iter = retrieve(params, meta)
    record_iter = align_vfrag(row_iter, params, meta)

    print(f"Try to insert data into {params['table_name']}...")
    bc = BuildComplex(params['verbose'])
    bc.insert_batch_records(record_iter, params['table_name'])

    footer(meta)
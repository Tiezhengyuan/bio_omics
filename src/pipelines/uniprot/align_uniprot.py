'''
requirements:
    - blastp alignment *.blastp5 files
    - table meta_blastp5
example: 
    - parallel -j8 bomics align_uniprot ::: {0..241}
    - bomics align_uniprot 0
functions:
    - pust alignment results of blastp into align_uniprot
'''

from src.helper import *
from bioomics import ParseBlastXml

@header
def retrieve_pdb(params, meta):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        SELECT B.pdb_id, B.relative_aln
        FROM meta_blastp5  B
        LEFT JOIN meta_pdb P ON B.pdb_id=P.pdb_id
        WHERE P.pdb_group = {pdb_group}
            AND B.pdb_id NOT IN (
                SELECT pdb_id FROM {table_name}
           )
    ;"""
    rows = QueryComplex(params['verbose']).list_data(query)
    meta['align_uniprot'] = len(rows)
    for row in rows:
        row['aln_file'] = params['output_pdb_dir'] / row['relative_aln']
        if row['aln_file'].is_file():
            yield row

def align_uniprot(row_iter, meta):
    for row in row_iter:
        best_hits = ParseBlastXml(row['aln_file']).best_hits()
        for hit in best_hits:
            rec = (row['pdb_id'], hit['query'])
            if 'hit_desc' in hit:
                uniprot_acc = hit['hit_desc'].split('|')[1]
                rec += (
                    uniprot_acc,
                    float(hit['evalue']),
                    int(hit['bit_score']),
                    float(hit['identity']),
                )
                meta['records'] += 1
            else:
                # insert null if not hits detected
                rec += (None,) * 4
                meta['null'] += 1
            yield rec
    
if __name__ == "__main__":
    params.update({
        'table_name': 'align_uniprot',
    })
    row_iter = retrieve_pdb(params, meta)
    record_iter = align_uniprot(row_iter, meta)

    print(f"Try to insert alignment into table {params['table_name']}")
    bc = BuildComplex(params['verbose'])
    bc.insert_batch_records(record_iter, params['table_name'])

    footer(meta)
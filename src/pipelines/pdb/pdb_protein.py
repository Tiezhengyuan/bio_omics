'''
requirements:
    - enough free space ~500GB
example:
    - parallel -j8 python app.py pdb_protein ::: {0..240}
    - python app.py pdb_protein 0
functions:
 - download fasta sequences
 - update table pdb_protein
 '''

from src.helper import *
from bioomics import PullPdb, ProcessFasta

@header
def retrieve(params, meta):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        SELECT pdb_id FROM meta_pdb
        WHERE pdb_group = {pdb_group}
        EXCEPT
        SELECT pdb_id FROM {table_name}
    """
    rows = QueryComplex(params['verbose']).list_data(query)
    meta['rows'] = len(rows)
    return rows

def download_fasta(rows, params, meta):
    p = PullPdb(params['rawdata_pdb_dir'], params['verbose'])
    for row in rows:
        pdb_id = row['pdb_id']
        fa_file, endpoint = p.download_fasta(pdb_id)
        if fa_file:
            records = ProcessFasta(fa_file, params['verbose']).load_fasta()
            for record in records:
                items = record.description.split('|')
                pro_id = items[0]
                s1 = re.findall(r'Chains? ([\w|,| ]*)\[?', items[1])
                chain = s1[0].replace(' ', '') if s1 else None
                s2 = re.findall(r'\[auth (.)\]', items[1])
                auth_chain = s2[0] if s2 else None
                pro_seq = str(record.seq)
                pro_len = len(pro_seq)
                meta['seqs'] += 1
                yield (pdb_id, endpoint, pro_id, chain, auth_chain, pro_len, pro_seq)
        else:
            meta['null'] += 1
            yield (pdb_id, ) + (None,) * 6


if __name__ == "__main__":
    params.update({
        'overwrite': False,
        'chunk_size': 10,
        'table_name': 'pdb_protein',
    })

    rows = retrieve(params, meta)
    record_iter = download_fasta(rows, params, meta)
    
    bc = BuildComplex(params['verbose'], params['chunk_size'])
    bc.insert_batch_records(record_iter, params['table_name'])

    footer(meta)
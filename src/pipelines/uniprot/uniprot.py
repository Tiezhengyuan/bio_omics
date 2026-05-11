'''
requirements:
    - uniprot_sprot.dat is downloaded
example:
    - python app.py uniprot
function:
    - synchronize PDB/UniProt/alignments with data in database COMPLEX.
'''

from src.helper import *
from bioomics import Uniprot, ParseDatRecord

def get_record(params, meta):
    qc = QueryComplex(params['verbose'])
    db_acc = qc.unique_values(params['table_name'], 'uniprot_acc')
    u = Uniprot(params['uniprot_sprot_dat'])
    for record in u.dat_iter():
        uniprot_acc = record.id
        if uniprot_acc not in db_acc:
            meta['uniprot_records'] += 1
            yield record
        else:
            meta['skip'] += 1
        
def build(record_iter, meta):
    for record in record_iter:
        p = ParseDatRecord(record)
        meta['insert_records'] += 1
        rec = (p.uniprot_acc, p.name, p.organism, len(p.seq), p.seq,)
        yield rec
    
##################################
if __name__ == "__main__":
    params.update({
        'uniprot_sprot_dat': os.getenv('uniprot_sprot_dat'),
        'table_name': 'uniprot',
    })
    # prepare record
    data_iter = get_record(params, meta)
    record_iter = build(data_iter, meta)

    print(f"Try to insert UniProt data into table {params['table_name']}...")
    bc = BuildComplex(params['verbose'])
    bc.insert_batch_records(record_iter, params['table_name'])

    footer(meta)


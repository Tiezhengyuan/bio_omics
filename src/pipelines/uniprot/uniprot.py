'''
requirements:
    - uniprot_sprot.dat is downloaded
example:
    - bomics uniprot
function:
    - synchronize PDB/UniProt/alignments with data in database COMPLEX.
'''

from src.helper import *
from bioomics import Uniprot, ParseDatRecord

def get_record(params, meta):
    qc = QueryComplex(params['verbose'])
    db_acc = qc.unique_values(params['table_name'], 'uniprot_acc')
    
    u = Uniprot(params['rawdata_uniprot_dir'])
    for record in u.filter_uniprot_sprot_dat(db_acc):
        meta['uniprot_records'] += 1
        yield record
        
def build(record_iter, meta):
    for record in record_iter:
        p = ParseDatRecord(record)
        meta['insert_records'] += 1
        rec = (
            p.uniprot_acc,
            p.name,
            p.organism,
            p.seq,
            p.creation_date,
            p.seq_date
        )
        yield rec
    
##################################
if __name__ == "__main__":
    params.update({
        'table_name': 'uniprot',
    })
    # prepare record
    data_iter = get_record(params, meta)
    record_iter = build(data_iter, meta)

    print(f"Try to insert UniProt data into table {params['table_name']}...")
    bc = BuildComplex(params['verbose'])
    bc.insert_batch_records(record_iter, params['table_name'])

    footer(meta)
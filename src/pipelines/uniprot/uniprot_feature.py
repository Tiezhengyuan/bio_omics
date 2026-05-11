'''
requirements:
    - uniprot_sprot.dat is downloaded
example:
    - python app.py uniprot_feature
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
            meta['uniprot_record'] += 1
            yield record
        else:
            meta['skip'] += 1
        
def build(data_iter, meta):
    for record in data_iter:
        tag = 0
        uniprot_acc = record.id
        features = ParseDatRecord(record).features()
        for feature in features:
            rec = (
                uniprot_acc,
                feature.get('type'),
                feature.get('start'),
                feature.get('end'),
                feature.get('seq'),
                feature.get('note'),
                feature.get('evidence'),
            )
            meta['insert_records'] += 1
            tag =1
            yield rec
        if tag == 0:
            meta['null'] += 1
            yield (uniprot_acc,) + (None,) * 6

##################################
if __name__ == "__main__":
    params.update({
        'uniprot_sprot_dat': os.getenv('uniprot_sprot_dat'),
        'table_name': 'uniprot_feature',
    })
    # prepare record
    data_iter = get_record(params, meta)
    record_iter = build(data_iter, meta)

    print(f"Try to insert UniProt data into table {params['table_name']}...")
    bc = BuildComplex(params['verbose'])
    bc.insert_batch_records(record_iter, params['table_name'])

    footer(meta)
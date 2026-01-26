'''
requirements:
    - uniprot_sprot.dat is downloaded
example:
    - python src/pipelines/uniprot_feature.py
function:
    - synchronize PDB/UniProt/alignments with data in database COMPLEX.
'''

from helper import *
from bioomics import Uniprot

def get_record(args, meta):
    qc = QueryComplex(args['verbose'])
    db_acc = qc.unique_values(args['table_name'], 'uniprot_acc')
    u = Uniprot(args['uniprot_sprot_dat'])
    for record in u.dat_iter():
        uniprot_acc = record.id
        if uniprot_acc not in db_acc:
            meta['record'] += 1
            yield record
        else:
            meta['skip'] += 1
        

##################################
if __name__ == "__main__":
    args = {
        'verbose': True,
        'chunk_size': 100,
        'uniprot_sprot_dat': os.getenv('uniprot_sprot_dat'),
        'table_name': 'uniprot_feature',
    }
    # prepare record
    records_iter = get_record(args, meta)

    print(f"Try to insert UniProt data into table {args['table_name']}...")
    bc = BuildComplex(records_iter, args['verbose'], args['chunk_size'])
    bc.insert_uniprot_feature()
    print(f"Synchronization data with {args['table_name']} is good.")



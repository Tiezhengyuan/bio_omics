'''
requirements:
    - uniprot_sprot.dat is downloaded
example:
    - python src/pipelines/put_uniprot_keyword.py
function:
    - synchronize PDB/UniProt/alignments with data in database COMPLEX.
'''
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import QueryComplex,  BuildComplex, Uniprot



##################################
if __name__ == "__main__":
    args = {
        'verbose': True,
        'overwrite': False,
        'uniprot_sprot_dat': '/home/yuan/data/UniProt/uniprot_sprot.dat.gz',
    }
    table_name = 'uniprot_keyword'
    print(f"Try to insert UniProt data into table {table_name}...")
    u = Uniprot(args['uniprot_sprot_dat'])

    # prepare iterator
    if args['overwrite']:
        records_iter = u.dat_iter()
        bc = BuildComplex(records_iter, args['verbose'])
        bc.empty_table(table_name)
        bc.insert_uniprot_keyword()
        print(f'Insert entire data into {table_name}.')
    else:
        qc = QueryComplex(args['verbose'])
        db_ids = qc.get_unique_values(table_name, 'uniprot_acc')
        records_iter = u.filter_dat_iter(db_ids)
        bc = BuildComplex(records_iter, args['verbose'])
        bc.insert_uniprot_keyword()
        print(f'Synchronization data with {table_name} is good.')



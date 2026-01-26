#!/usr/bin/python
'''
requirements:
    - PDB raw data should be downloaded to the local
example:
    - python src/pipelines/pdb.py
function:
    - build table pdb and insert data into table pdb
'''

from helper import *
from bioomics import PullPdb

def scan_raw_data(db_pdb_ids, args, meta):
    raw_pdb_dir = os.path.join(args['pdb_dir'], args['data_format'])
    pdb_path_dict = PullPdb(raw_pdb_dir).get_raw_pdb()

    ud = UpdateComplex(args['table_name'], args['verbose'])
    for pdb_id in db_pdb_ids:
        if pdb_id in pdb_path_dict:
            del pdb_path_dict[pdb_id]
        else:
            res = ud.delete_rows('pdb_id', pdb_id)
            if res:
                meta['delete_pdb'] += 1
            else:
                meta['delete_fail'] += 1
    meta['new_pdb_files'] = len(pdb_path_dict)
    return pdb_path_dict


if __name__ == "__main__":
    args = {
        'verbose': True,
        'pdb_dir': os.getenv('pdb_dir'),
        'data_format': 'pdb',
        'table_name': 'pdb',
    }
    # get pdb in db
    qc = QueryComplex(args['verbose'])
    db_pdb_ids = qc.unique_values(args['table_name'], 'pdb_id')
    meta['db_pdb'] = len(db_pdb_ids)

    # scan raw data
    pdb_path_dict = scan_raw_data(db_pdb_ids, args, meta)
    pdb_files = pdb_path_dict.values()
    BuildComplex(pdb_files, args['verbose']).insert_pdb()

    print(meta)
    print("Done. Great\n\n")
#!/usr/bin/python
'''
requirements:
    - PDB raw data should be downloaded to the local
example:
    - python app.py pdb
function:
    - build table pdb and insert data
'''

from src.helper import *
from bioomics import PullPdb

def scan_raw_data(db_pdb_ids, params, meta):
    raw_pdb_dir = os.path.join(params['pdb_dir'], params['data_format'])
    pdb_path_dict = PullPdb(raw_pdb_dir).get_raw_pdb()

    ud = UpdateComplex(params['table_name'], params['verbose'])
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
    params = {
        'verbose': True,
        'pdb_dir': os.getenv('pdb_dir'),
        'data_format': 'pdb',
        'table_name': 'pdb',
    }
    # get pdb in db
    qc = QueryComplex(params['verbose'])
    db_pdb_ids = qc.unique_values(params['table_name'], 'pdb_id')
    meta['db_pdb'] = len(db_pdb_ids)

    # scan raw data
    pdb_path_dict = scan_raw_data(db_pdb_ids, params, meta)
    pdb_files = pdb_path_dict.values()
    BuildComplex(pdb_files, params['verbose']).insert_pdb()

    footer(meta)

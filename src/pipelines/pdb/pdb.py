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
from bioomics import PullPdb, Pdb

def scan(db_pdb_ids, params, meta):
    pdb_path_dict = PullPdb(params['rawdata_pdb_dir']).get_raw_pdb()

    dc = DeleteComplex(params['verbose'])
    for pdb_id in db_pdb_ids:
        if pdb_id in pdb_path_dict:
            del pdb_path_dict[pdb_id]
        else:
            res = dc.delete_rows(params['table_name'], 'pdb_id', [pdb_id,])
            if res:
                meta['delete_pdb'] += 1
            else:
                meta['delete_fail'] += 1
    meta['new_pdb_files'] = len(pdb_path_dict)
    pdb_files = pdb_path_dict.values()
    return pdb_files

def build(pdb_files, params, meta):
    for pdb_file in pdb_files:
        p = Pdb(pdb_file, None, params['verbose'])
        p.load_structure()

        pdb_id = p.structure_id
        pdb_url = f"https://www.rcsb.org/structure/{pdb_id}"
        release_date = p.structure.header['release_date']
        name = p.structure.header['name']
        rec = (pdb_id, pdb_url, release_date, name)
        meta['records'] += 1
        yield rec

if __name__ == "__main__":
    params.update({
        'table_name': 'pdb',
    })
    # get pdb in db
    qc = QueryComplex(params['verbose'])
    db_pdb_ids = qc.unique_values(params['table_name'], 'pdb_id')
    meta['db_pdb'] = len(db_pdb_ids)

    # scan raw data
    pdb_files = scan(db_pdb_ids, params, meta)
    record_iter = build(pdb_files, params, meta)

    # update table pdb
    bc = BuildComplex(params['verbose'])
    bc.insert_batch_records(record_iter, params['table_name'])

    footer(meta)

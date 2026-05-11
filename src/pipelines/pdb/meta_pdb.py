#!/usr/bin/python
'''
requirements:
    - PDB raw data should be downloaded to the local
    - table pdb
example:
    - python app.py meta_pdb
function:
    - build table pdb and insert data into table meta_pdb
'''
from src.helper import *
from bioomics import ProcessJson

def get_meta_json(pdb_ids, params, meta):
    list_dir = os.path.join(params['output_pdb_dir'], '_list', params['data_format'])
    for file_name in os.listdir(list_dir):
        if file_name.endswith('.json'):
            pdb_group = os.path.splitext(file_name)[0]
            path = os.path.join(list_dir, file_name)
            print(f"Retrieve meta data from {path}")
            pdb_meta = ProcessJson.load_json(path)
            meta['pdb_group'] += 1
            for item in pdb_meta:
                pdb_id = item['pdb_id']
                if pdb_id not in pdb_ids:
                    relative_path = item['pdb_file'].replace(params['rawdata_pdb_dir'], '')
                    print(relative_path)
                    if relative_path.startswith('/'):
                        relative_path = relative_path[1:]
                    meta['new_pdb'] += 1
                    yield (pdb_id, relative_path, int(pdb_group))
                else:
                    meta['skip_pdb'] += 1

if __name__ == "__main__":
    params.update({
        'data_format': 'pdb',
        'table_name': 'meta_pdb',
    })
    # retrieve pdb in database
    qc = QueryComplex(params['verbose'])
    pdb_ids = qc.unique_values(params['table_name'], 'pdb_id')
    record_iter = get_meta_json(pdb_ids, params, meta)

    print(f"Try to insert PDB meta data into table {params['table_name']}...")
    bc = BuildComplex(params['verbose'])
    bc.insert_batch_records(record_iter, params['table_name'])
        
    footer(params, meta)

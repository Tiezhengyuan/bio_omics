#!/usr/bin/python
'''
requirements:
    - PDB raw data should be downloaded to the local
    - table pdb
example:
    - python src/pipelines/meta_pdb.py
function:
    - build table pdb and insert data into table meta_pdb
'''
from helper import *
from bioomics import ProcessJson

def get_meta_json(pdb_ids, args, meta):
    list_dir = os.path.join(args['data_dir'], '_list', args['data_format'])
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
                    relative_path = item['pdb_file'].replace(args['pdb_dir'], '')
                    if relative_path.startswith('/'):
                        relative_path = relative_path[1:]
                    meta['new_pdb'] += 1
                    yield (pdb_id, relative_path, int(pdb_group))
                else:
                    meta['skip_pdb'] += 1

if __name__ == "__main__":
    args = {
        'verbose': True,
        'pdb_dir': os.getenv("pdb_dir"),
        'data_dir': os.getenv("data_dir"),
        'data_format': 'pdb',
        'table_name': 'meta_pdb',
    }
    # retrieve pdb in database
    pdb_ids = QueryComplex(args['verbose']).unique_values(args['table_name'], 'pdb_id')
    record_iter = get_meta_json(pdb_ids, args, meta)

    print(f"Try to insert PDB meta data into table {args['table_name']}...")
    bc = BuildComplex(record_iter, args['verbose'])
    bc.insert_batch_records(args['table_name'])
        
    print(meta)
    print("Done. Great\n\n")

#!/usr/bin/python
'''
requirements:
    - PDB raw data should be downloaded to the local
    - table pdb
example:
    - python src/pipelines/put_meta_pdb.py
function:
    - build table pdb and insert data into table meta_pdb
'''
from helper import *
from bioomics import ProcessJson


if __name__ == "__main__":
    args = {
        'verbose': True,
        'pdb_dir': os.getenv("pdb_dir"),
        'data_dir': os.getenv("data_dir"),
        'data_format': 'pdb',
        'table_name': 'meta_pdb',
    }
    # retrieve pdb in database
    table_name = args['table_name']
    pdb_ids = QueryComplex(args['verbose']).unique_values(table_name, 'pdb_id')

    list_dir = os.path.join(args['data_dir'], '_list', args['data_format'])
    for file_name in os.listdir(list_dir):
        if file_name.endswith('.json'):
            pdb_group = os.path.splitext(file_name)[0]
            path = os.path.join(list_dir, file_name)
            print(f"Retrieve meta data from {path}")
            meta = ProcessJson.load_json(path)
            records = []
            for item in meta:
                pdb_id = item['pdb_id']
                if pdb_id not in pdb_ids:
                    relative_path = item['pdb_file'].replace(args['pdb_dir'], '')
                    if relative_path.startswith('/'):
                        relative_path = relative_path[1:]
                    rec = (pdb_id, relative_path, int(pdb_group))
                    records.append(rec)
            # insertion
            if records:
                print(f"Try to insert {len(records)} data into table {table_name}...")
                bc = BuildComplex(records, args['verbose'])
                bc.insert_batch_records(table_name)
    
    print("Done. Great\n\n")

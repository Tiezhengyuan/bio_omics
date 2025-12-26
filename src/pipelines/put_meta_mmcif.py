#!/usr/bin/python
'''
requirements:
    - PDB raw data should be downloaded to the local
    - table pdb
example:
    - python src/pipelines/put_meta_mmcif.py
function:
    - build table pdb and insert data into table meta_mmcif
'''
from dotenv import load_dotenv
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import ProcessJson, BuildComplex, QueryComplex


if __name__ == "__main__":
    load_dotenv()
    args = {
        'verbose': True,
        'pdb_dir': os.getenv("pdb_dir"),
        'data_dir': os.getenv("data_dir"),
        'data_format': 'mmCIF',
        'table_name': 'meta_mmcif',
    }
    
    list_dir = os.path.join(args['data_dir'], '_list', args['data_format'])
    # build reords
    table_name = args['table_name']
    pdb_ids = QueryComplex(args['verbose']).unique_values(table_name, 'pdb_id')
    for file_name in os.listdir(list_dir):
        if file_name.endswith('.json'):
            pdb_group = os.path.splitext(file_name)[0]
            path = os.path.join(list_dir, file_name)
            print(f"Load meta data from {path}")
            meta = ProcessJson.load_json(path)
            records = []
            for item in meta:
                pdb_id = item['pdb_id']
                print
                if pdb_id not in pdb_ids:
                    relative_path = item['pdb_file'].replace(args['pdb_dir'], '')
                    rec = (pdb_id, relative_path, int(pdb_group))
                    records.append(rec)
            # insertion
            if records:
                print(f"Try to insert {len(records)} data into table {table_name}...")
                bc = BuildComplex(records, args['verbose'])
                bc.insert_batch_records(table_name)
        break
    
    print("Done. Great\n\n")

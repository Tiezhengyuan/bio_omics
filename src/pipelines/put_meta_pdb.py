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
import os
import re
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import ProcessJson, BuildComplex, QueryComplex


if __name__ == "__main__":
    args = {
        'verbose': True,
        'list_dir': '/home/yuan/data/pdb/_list'

    }

    # build reords
    table_name = 'meta_pdb'
    pdb_ids = QueryComplex().unique_values(table_name, 'pdb_id')
    for file_name in os.listdir(args['list_dir']):
        if file_name.endswith('.json'):
            pdb_group = re.split(r'_|\.', file_name)[1]
            path = os.path.join(args['list_dir'], file_name)
            print(f"Load meta data from {path}")
            meta = ProcessJson.load_json(path)
            records = []
            for item in meta:
                pdb_id = item['pdb_id']
                if pdb_id not in pdb_ids:
                    pdb_file = item['pdb_file']
                    rec = (pdb_id, pdb_file, int(pdb_group))
                    records.append(rec)
            # insertion
            if records:
                print(f"Try to insert {len(records)} data into table {table_name}...")
                bc = BuildComplex(records, args['verbose'])
                bc.insert_meta_pdb()
    
    print("Done. Great\n\n")

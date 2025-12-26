'''
Download data from database SAbDab https://opig.stats.ox.ac.uk
example: 
    - python src/pipelines/ab_sabdab.py
functions:
    - download data from SAbDab
    - put data to table ab_sabdab
'''
from collections import defaultdict
from dotenv import load_dotenv
import numpy as np
import pandas as pd
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import BuildComplex, QueryComplex

def prepare_data(rows):
    for row in rows:
        pdb_id = row['pdb_id']
        path = row['summary_file']
        df = pd.read_csv(path, sep='\t')
        for _, row in df.iterrows():
            row = row.replace({
                np.nan: None,
                'NOT': None,
            })
            if row.get('pdb'):
                pdb_id = str(row['pdb']).upper()
                resolution = row.get('resolution')
                if resolution is None:
                    resolution = 0
                elif isinstance(resolution, str):
                    resolution = float(resolution.split(',')[0])
                record = (
                    pdb_id,
                    row['model'],
                    row['Hchain'],
                    row['Lchain'],
                    row['antigen_chain'],
                    row['antigen_type'],
                    resolution,
                )
                yield record
        

if __name__ == '__main__':
    meta = defaultdict(int)
    load_dotenv()
    args = {
        'verbose': True,
        'chunk_size': 5,
        'download_dir': os.getenv("download_dir"),
        'table_name': 'ab_sabdab',
    }
    # 
    query = "SELECT * FROM meta_sabdab WHERE summary_file IS NOT NULL;"
    rows = QueryComplex(args['verbose']).list_data(query)
    meta['rows'] = len(rows)

    # insertion
    
    record_iter = prepare_data(rows)
    bc = BuildComplex(record_iter, args['verbose'], args['chunk_size'])
    bc.empty_table(args['table_name'])
    bc.insert_batch_records(args['table_name'])

    print(meta)

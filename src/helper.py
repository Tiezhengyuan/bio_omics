'''
set imports required by all pipelines
'''
from collections import defaultdict
from datetime import datetime
from dotenv import load_dotenv
import functools
import itertools
import os
import re
from pathlib import Path
import sys
import numpy as np
import pandas as pd

src_dir = os.path.dirname(__file__)
if src_dir not in sys.path:
    sys.path.append(src_dir)

env_file = os.path.join(src_dir, '.env')
load_dotenv(env_file)

meta = defaultdict(int)

import socket
from bioomics.common import Dir, WriteTmp, ProcessPickle
from bioomics.database import BuildComplex, QueryComplex, UpdateComplex
from bioomics.PDB import Pdb
from bioomics.protein import ProcessFasta

def header(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        params, meta = args
        pdb_group = params.get('pdb_group')
        table = params.get('table_name')
        print("#########\n")
        print(f"pdb group = {pdb_group}")
        print(f"Try to retrieve data from {table}...")
        meta['pdb_group'] = pdb_group
        meta['host'] = socket.gethostname()
        result = func(*args, **kwargs)
        return result
    return wrapper

def validate_data(col_name, params, meta):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    # scan no existing files
    print(f"Scan {table_name}.{col_name}, group={pdb_group}.", end=' ')
    query = f"""
        SELECT S.pdb_id, S.{col_name}
        FROM {table_name}  S
        LEFT JOIN meta_pdb P ON S.pdb_id = P.pdb_id
        WHERE P.pdb_group = {pdb_group}
            AND S.{col_name} IS NOT NULL
    ;"""
    rows = QueryComplex(params['verbose']).list_data(query)
    del_pool = []
    for row in rows:
        pdb_file = os.path.join(params['data_dir'], row[col_name])
        if not os.path.isfile(pdb_file):
            del_pool.append(row[col_name])
    
    # delete rows
    n = len(del_pool)
    if n > 0:
        print(f"Number of {n} out of {len(rows)} are not existing.", end=' ')
        print(f"Try to delete rows in {table_name}...", end='')
        uc = UpdateComplex(table_name, params['verbose'])
        chunk_size = 100
        for i in range(0, n-chunk_size, chunk_size):
            values = del_pool[i:i+chunk_size]
            uc.delete_rows(col_name, values)
            meta['delete_rows'] += len(values)
            print(i, end=' ')
        print(f'{n}\n')

def footer(meta):
    print('Statistics of this run', meta)
    print("Done!\n\n")
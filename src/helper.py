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
import shutil
import sys
import numpy as np
import pandas as pd
import subprocess
import socket

src_dir = os.path.dirname(__file__)
if src_dir not in sys.path:
    sys.path.append(src_dir)

env_file = os.path.join(src_dir, '.env')
load_dotenv(env_file)

# initiate parameters
params = {
    'pdb_group': sys.argv[1] if len(sys.argv) > 1 else None,
    'verbose': True if os.getenv('verbose') == 'true' else False,
    'pdb_dir': os.getenv('pdb_dir'),
    'data_dir': os.getenv('data_dir'),
    'simulate_dir': os.getenv('simulate_dir'),
    'predict_dir': os.getenv('predict_dir'),
}
meta = defaultdict(int)

import socket
from bioomics.common import Dir, WriteTmp, ProcessPickle, ProcessJson
from bioomics.database import BuildComplex, UpdateComplex, QueryComplex
from bioomics.PDB import Pdb
from bioomics.protein import ProcessFasta


def local_ip():
    s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    s.connect(("8.8.8.8", 80))
    ip_addr = s.getsockname()[0]
    s.close()
    return ip_addr

def add_params(params:dict={}):
    '''
    update args and meta
    '''
    # host name
    params['host'] = socket.gethostname()
    # IP address of this host
    params['ip_addr'] = local_ip()
    params['script'] = sys.argv[0]

    print('Try to assign pdb_group...')
    if len(sys.argv) > 1:
        # by arguments in command line
        params['pdb_group'] = sys.argv[1]
    else:
        #  by table meta_execute
        pdb_group = QueryComplex().pop_pdb_group(params['table_name'])
        params['pdb_group'] = pdb_group
        if params['pdb_group'] is not None:
            # add new record
            rec = (params['table_name'], params['pdb_group'], params['host'], params['ip_addr'], params['script'])
            BuildComplex([rec,], params['verbose']).insert_batch_records('meta_execute')
            params['meta_execute'] = True
            if params['verbose']:
                print('Update meta_exeucte:', rec)

    if params['pdb_group'] is None:
        print(f"\nExit! Try process {params['table_name']}, but not pdb_group is assigned.")
        sys.exit(0)
    if params['verbose']:
        print('Parameters', params)
    return params

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

def copy_data(outdir):
    data_hostname = os.getenv('DATA_HOSTNAME')
    execute_hostname = socket.gethostname()
    if data_hostname == execute_hostname:
        print(f"Skip synchronization, identical IP={data_hostname}")
    else:
        print(f"Try to copy {outdir} from {execute_hostname} to {data_hostname}.")
        data_host_user = os.getenv('DATA_HOST_USER')
        data_host_ip = os.getenv('DATA_HOST_IP')
        cmd = ['rsync', '-avz', '--progress',
            '--mkpath', '--remove-source-files',
            f"{outdir}/*",
            f"{data_host_user}@{data_hostname}:{outdir}/",
        ]
        try:
            cmd = ' '.join(cmd)
            print(cmd)
            res = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
            print('Synchronize results: ', res.stdout, res.stderr)
        except Exception as e:
            print(f"ERROR: Synchronization between {outdir}, err={e}")

def footer(meta, params:dict={}):
    # update table meta_execute with end time
    if params.get('meta_execute'):
        UpdateComplex('meta_execute', params['verbose']).execute_end_time(params)
    # print statistics of this run
    print('Statistics of this run', meta)
    print("Done!\n\n")
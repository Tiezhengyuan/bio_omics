'''
set imports required by all pipelines
'''
from collections import defaultdict
from datetime import datetime
from dotenv import load_dotenv
import functools
import itertools
import numpy as np
import os
import pandas as pd
from pathlib import Path
import re
import shutil
import socket
import sys
import subprocess
import torch

# load env
src_dir = os.path.dirname(__file__)
if src_dir not in sys.path:
    sys.path.append(src_dir)

project_dir = os.path.dirname(src_dir)
env_file = os.path.join(project_dir, '.env')
load_dotenv(env_file)

# load customary modules
from bioomics.database import BuildComplex, UpdateComplex, QueryComplex, DeleteComplex


# initiate parameters
params = {
    'host': socket.gethostname(),
    'pdb_group': sys.argv[1] if len(sys.argv) > 1 else None,
    'verbose': True if os.getenv('verbose') == 'true' else False,
    'rawdata_dir': Path(os.getenv('rawdata_dir')),
    'output_dir': Path(os.getenv('output_dir')),
    'tools_dir': Path(os.getenv('tools_dir')),
    'gpu': torch.cuda.is_available(),
}
params['rawdata_pdb_dir'] = params['rawdata_dir'] / 'pdb'
params['rawdata_uniprot_dir'] = params['rawdata_dir'] / 'UniProt'
params['output_pdb_dir'] = params['output_dir'] / 'pdb'
params['output_protein_dir'] = params['output_dir'] / 'protein'
params['simulate_dir'] = params['output_dir'] / 'simulate'
params['predict_dir'] = params['output_dir'] / 'predict'
params['design_dir'] = params['output_dir'] / 'design'
meta = defaultdict(int)


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
            bc = BuildComplex(params['verbose'])
            bc.insert_batch_records([rec,], 'meta_execute')
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



#TODO: deprecated in the future
def validate_data(col_name, params, meta, outdir:Path=None):
    if outdir is None:
        outdir = params['output_pdb_dir']
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    # scan no existing files
    print(f"Scan {table_name}.{col_name}, group={pdb_group}...")
    query = f"""
        SELECT S.pdb_id, S.{col_name}
        FROM {table_name}  S
        LEFT JOIN meta_pdb P ON S.pdb_id = P.pdb_id
        WHERE P.pdb_group = {pdb_group} AND S.{col_name} IS NOT NULL
    ;"""
    rows = QueryComplex(params['verbose']).list_data(query)
    del_pool = []
    for row in rows:
        pdb_file = outdir / row[col_name]
        if not pdb_file.is_file():
            print(pdb_file)
            del_pool.append(row[col_name])
    
    # delete rows
    n = len(del_pool)
    if n > 0:
        print(f"Number of {n} out of {len(rows)} are not existing.", end=' ')
        dc = DeleteComplex(params['verbose'])
        chunk_size = 100
        end = n - chunk_size if n > chunk_size else n
        for i in range(0, end, chunk_size):
            values = del_pool[i:i+chunk_size]
            res = dc.delete_rows(table_name, col_name, values)
            if res:
                meta['delete_rows'] += len(values)
                print(i, end=' ')
        print(f'{n}\n')
    else:
        print(f"PASS: Validate {table_name}.{col_name}.")

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
        uc = UpdateComplex(params['verbose'])
        uc.execute_end_time(params)
    # print statistics of this run
    print('Statistics of this run', meta)
    print("Done!\n\n")
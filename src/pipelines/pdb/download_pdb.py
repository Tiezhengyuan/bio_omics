'''
requirements:
    - enough free space ~500GB
example:
    - python app.py download_pdb
functions:
 - download the latest data in pdb and mmcif format
 - wrap the commands for running pipelines in multi-rocessing
'''

from src.helper import *
from bioomics import ProcessJson, PullPdb

if __name__ == "__main__":
    params.update({
        'overwrite': False,
        'data_format': 'pdb',
    })
    
    # donwload pdb
    p = PullPdb(params['rawdata_dir'], params['verbose'])
    p.synchronize_pdb(params['data_format'])

    # prepare _list dir
    list_dir = params['output_pdb_dir'] /'_list' / params['data_format']
    print(list_dir)
    shutil.rmtree(list_dir, ignore_errors=True)
    list_dir.mkdir(parents=True, exist_ok=True)
    
    # build json
    raw_dir = params['rawdata_dir'] / params['data_format']
    pdb_iter = PullPdb(raw_dir, params['verbose']).iter_raw_pdb_args(params)
    file_list = ProcessJson.iter_to_json(pdb_iter, list_dir, chunk_size=1_000)
    print(f"\n{file_list}\n")

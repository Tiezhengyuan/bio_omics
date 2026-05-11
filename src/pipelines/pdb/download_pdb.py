'''
requirements:
    - enough free space ~500GB
example:
    - python app.py download_pdb
functions:
 - download the latest pdb files
 - wrap the commands for running pipelines in multi-rocessing
'''

from src.helper import *
from bioomics import Dir, ProcessJson, PullPdb

if __name__ == "__main__":
    params.update({
        'overwrite': False,
        'data_format': 'pdb',
    })
    
    # donwload pdb
    raw_dir = os.path.join(params['rawdata_pdb_dir'], params['data_format'])
    PullPdb(raw_dir, params['verbose']).synchronize_pdb(params['data_format'])

    # prepare _list dir
    list_dir = os.path.join(params['output_pdb_dir'], '_list', params['data_format'])
    d = Dir(list_dir)
    d.init_dir()
    d.clear_dir()
    
    # build json
    pdb_iter = PullPdb(raw_dir).iter_raw_pdb_args(params)
    file_list = ProcessJson.iter_to_json(pdb_iter, list_dir, chunk_size=1_000)
    print(f"\n{file_list}\n")

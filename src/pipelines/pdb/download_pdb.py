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
from bioomics import ProcessJson, PullPdb

if __name__ == "__main__":
    args = {
        'verbose': True,
        'overwrite': False,
        'pdb_dir': os.getenv("pdb_dir"),
        'data_dir': os.getenv("data_dir"),
        'data_format': 'pdb',
    }
    # donwload pdb
    raw_dir = os.path.join(args['pdb_dir'], args['data_format'])
    PullPdb(raw_dir, args['verbose']).synchronize_pdb(args['data_format'])

    # prepare _list dir
    list_dir = os.path.join(args['data_dir'], '_list', args['data_format'])
    d = Dir(list_dir)
    d.init_dir()
    d.clear_dir()
    
    # build json
    pdb_iter = PullPdb(raw_dir).iter_raw_pdb_args(args)
    file_list = ProcessJson.iter_to_json(pdb_iter, list_dir, chunk_size=1_000)
    print(f"\n{file_list}\n")

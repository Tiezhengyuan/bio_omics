'''
requirements:
    - enough free space ~500GB
example:
    - python src/pipelines/download_mmcif.py
functions:
 - download the latest pdb files
 - wrap the commands for running pipelines in multi-rocessing
'''
from dotenv import load_dotenv
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import PullPdb, ProcessJson, Dir


if __name__ == "__main__":
    load_dotenv()
    args = {
        'verbose': True,
        'overwrite': False,
        'pdb_dir': os.getenv("pdb_dir"),
        'data_dir': os.getenv("data_dir"),
        'data_format': 'mmCIF',
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
    file_list = ProcessJson.iter_to_json(pdb_iter, list_dir, chunk_size=2_000)
    print(f"\n{file_list}\n")

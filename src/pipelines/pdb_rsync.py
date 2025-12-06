'''
requirements:
    - enough free space ~500GB
example:
    - python src/pipelines/pdb_rsync.py
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
    }
    # donwload pdb
    raw_dir = os.path.join(args['pdb_dir'], 'rsync')
    PullPdb(raw_dir, args['verbose']).synchronize_pdb()

    # prepare _list dir
    list_dir = os.path.join(args['data_dir'], '_list')
    d = Dir(list_dir)
    d.init_dir()
    d.clear_dir()
    outprefix = os.path.join(list_dir, 'pdb')
    
    # build json
    pdb_iter = PullPdb(raw_dir).iter_raw_pdb_args(args)
    file_list = ProcessJson.iter_to_json(pdb_iter, outprefix, chunk_size=2_000)
    print(f"\n{file_list}\n")

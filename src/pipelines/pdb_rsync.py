'''
example: python src/pipelines/pdb_rsync.py
functions:
 - download the latest pdb files
 - wrap the commands for running pipelines in multi-rocessing
'''
from Bio import SeqIO
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import PullPdb, ProcessJson, Dir


def wrap_args(args):
    # prepare _list dir
    args_dir = os.path.join(args['data_dir'], '_list')
    d = Dir(args_dir)
    d.init_dir()
    d.clear_dir()

    # build json
    pdb_iter = PullPdb(args['pdb_dir']).iter_raw_pdb_args(args)
    outprefix = os.path.join(args_dir, 'pdb')
    file_list = ProcessJson.iter_to_json(pdb_iter, outprefix, chunk_size=5_000)
    print('\n\n')
    for script in ('pdb_outfile.py', 'pdb_blast.py', 'pdb_igblastp.py',):
        print(f"parallel -j6 python src/pipelines/{script} ::: {args_dir}/pdb_*.json")
    print('\n\n')

if __name__ == "__main__":
    args = {
        'overwrite': False,
        'pdb_dir': '/home/yuan/data/pdb/rsync',
        'data_dir': '/home/yuan/data/pdb',
    }
    # donwload pdb
    PullPdb(args['pdb_dir']).synchronize_pdb()
    # wrap the commands for running pipelines
    wrap_args(args)
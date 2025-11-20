'''
database SAbDab https://opig.stats.ox.ac.uk
'''
from collections import defaultdict
import os
import numpy as np
import pandas as pd
import requests
import subprocess

from ..common import Dir

class Sabdab:
    url = 'https://opig.stats.ox.ac.uk/webapps/abdb/entries'

    def __init__(self, download_dir:str, verbose:bool=False, overwrite:bool=False):
        self.outdir = os.path.join(download_dir, 'SAbDab')
        Dir(self.outdir).init_dir()
        self.summary_dir = os.path.join(self.outdir, 'summary')
        Dir(self.summary_dir).init_dir()
        self.verbose = verbose
        self.overwrite = overwrite

    def retrieve_sabdab_pdb_ids(self):
        '''
        retrieve list of pdb_id
        '''
        try:
            # https://opig.stats.ox.ac.uk/webapps/abdb/entries/sabdab_pdbs_list.txt
            res = requests.get(f'{self.url}/sabdab_pdbs_list.txt')
            pdb_ids = res.text.split('\n')
            return pdb_ids
        except Exception as e:
            print(f"Error: can't retrieve pdb ids from {self.url}. error={e}")
        return []

    def download_summary(self, pdb_ids:list, meta:defaultdict=None):
        if meta is None:
            meta = defaultdict(int)
        if not pdb_ids:
            print(f"WARNING: summary can not be downloaded.")
            return None

        with open(os.path.join(self.outdir, 'error_summary.txt'), 'w') as f:
            for pdb_id in pdb_ids:
                outfile = os.path.join(self.summary_dir, f'{pdb_id}.tsv')
                if self.overwrite or not os.path.isfile(outfile):
                    endpoint = f'{self.url}/{pdb_id}/summary/{pdb_id}.tsv'
                    cmd = ['wget', '-c', endpoint, '-P', self.summary_dir]
                    try:
                        subprocess.run(cmd, check=True)
                        meta['succeed'] += 1
                    except Exception as e:
                        err = f"{pdb_id} | {e}\n"
                        f.write(err)
                        meta['fail'] += 1
                else:
                    meta['skip'] += 1
                pdb_id = str(pdb_id).upper()
                if pdb_id:
                    rec = (pdb_id, outfile) if os.path.isfile(outfile) else (pdb_id, None)
                    yield rec

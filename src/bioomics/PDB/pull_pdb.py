'''
process large number of PDB files 
'''
import json
import os
import re
import requests
import urllib
from datetime import date


from ..common.dir import Dir

class PullPdb:
    def __init__(self, pdb_dir:str=None, verbose=False):
        self.pdb_dir = pdb_dir
        self.verbose = verbose

    def synchronize_pdb(self, data_format:str=None):
        '''
        download entire PDB files to the local
        file_format: pdb, mmCIF, xml
        '''
        data_format = 'pdb' if data_format is None else data_format
        url = "rsync.rcsb.org"
        endpoint = f"ftp_data/structures/divided/{data_format}/"
        # -r recursive, -l copy symlinks as symlinks
        # -p keep permissions, -t keep time modifications
        # -v verbose, -z compress file
        cmd = f"rsync -rlpt -v -z --port=33444 {url}::{endpoint} {self.pdb_dir}"
        if self.verbose:
            print('Synchronize PDB data:', cmd)
        res = os.system(cmd)
        return res

    def download_pdb(self, pdb_id:str) -> str:
        '''
        download a certain pdb from PDB
        the local dir is ../pdb in default
        '''
        pdb_file = os.path.join(self.pdb_dir, f'{pdb_id}.pdb')
        if not os.path.isfile(pdb_file):
            print(f"download {pdb_id} to {pdb_file}")
            url = 'https://files.rcsb.org/view'
            os.system(f"wget -qnc {url}/{pdb_id}.pdb -P {pdb_dir}")
        print(f"{pdb_file} is ready.")
        return pdb_file

    def getpdb(self, pdb_entry):
        """
        Get the PDB file from SAbDab.
        Check that it has successfully downloaded.
        """
        out_file = os.path.join(self.pdb_dir, f"{pdb_entry}.pdb")
        endpoint = "https://opig.stats.ox.ac.uk/webapps/abdb/entries"
        url = f"{endpoint}/{pdb_entry}/structure/{pdb_entry}.pdb"
        urllib.request.urlretrieve(url, out_file)
        if os.path.isfile(out_file):
            Retrieved = open(out_file).read()
            if not Retrieved.count("ATOM"):
                print("Failed to retrieve PDB file from SAbDab")
                os.remove(out_file)
                return False
            else:
                return True
        else:
            return False

    def pdb_to_uniprot(self, pdb_id) -> dict:
        '''
        Uniprot-SIFTS
        mapping of pdb -> uniprot
        '''
        pdb_id = pdb_id.lower()
        url = f"https://www.ebi.ac.uk/pdbe/api/mappings/{pdb_id}"
        if self.verbose:
            print(f'Try to download pdb {pdb_id} from {url}.')
        try:
            res = requests.get(url)
            data = json.loads(res.text)
            return data[pdb_id]
        except Exception as e:
            print(e)
        return None

    def pdb_id(self, path):
        file_name = os.path.basename(path)
        if path.endswith('.ent.gz'):
            pdb_id = re.sub(r'^pdb|\.ent\.gz$', '', file_name)
            return pdb_id.upper()
        elif path.endswith('.cif.gz'):
            pdb_id = file_name.split('.')[0]
            return pdb_id.upper()
        elif path.endswith('.pdb'):
            pdb_id = re.sub(r'\pdb$', '', file_name)
            return pdb_id.upper()
        return file_name
    
    def iter_raw_pdb_args(self, args:dict):
        '''
        scan pdf files incursively
        '''
        file_iter = Dir(self.pdb_dir).recursive_files()
        for path in file_iter:
            pdb_id = self.pdb_id(path)
            if pdb_id:
                rec_args = dict(args)
                rec_args.update({
                    'pdb_id': pdb_id,
                    'pdb_file': path,
                })
                yield rec_args

    def iter_raw_pdb(self):
        '''
        scan pdf files incursively
        '''
        file_iter = Dir(self.pdb_dir).recursive_files()
        for path in file_iter:
            pdb_id = self.pdb_id(path)
            if pdb_id:
                yield pdb_id, path

    def get_raw_pdb(self) -> dict:
        res = {}
        file_iter = self.iter_raw_pdb()
        for pdb_id, path in file_iter:
            res[pdb_id] = path
        return res

    def filter_out_raw_pdb(self, pdb_ids) -> dict:
        n, res = 0, {}
        file_iter = self.iter_raw_pdb()
        for pdb_id, path in file_iter:
            n += 1
            if pdb_id not in pdb_ids: 
                res[pdb_id] = path
        print(f"Scan {n} PDB files in {self.pdb_dir}. Detect unknown {len(res)} PDB.")
        return res

    def filter_in_raw_pdb(pdb_ids):
        file_iter = self.iter_raw_pdb()
        for pdb_id, path in file_iter:
            if pdb_id in pdb_ids:
                yield pdb_id, path

    def get_new_raw_pdb(self, date_stamp) -> dict:
        res = {}
        file_iter = self.iter_raw_pdb()
        for pdb_id, path in file_iter:
            mod_time = os.path.getmtime(path)
            local_time = date.fromtimestamp(mod_time)
            if local_time >= date_stamp:
                res[pdb_id] = path
        return res

    def iter_meta(self):
        file_iter = Dir(self.pdb_dir).recursive_files()
        for path in file_iter:
            if path.endswith('meta.json'):
                pdb_id = os.path.basename(os.path.dirname(path))
                with open(path, 'r') as f:
                    data = json.load(f)
                    yield pdb_id, data


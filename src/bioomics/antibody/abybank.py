'''
http://www.abybank.org/
'''
from bs4 import BeautifulSoup
import re
import requests


class Abybank:
    url = 'http://www.abybank.org/abdb/Data'

    def __init__(self, verbose:bool=False):
        self.verbose = verbose
    
    def pull_pdb_records(self):
        abYbank_pdb = []
        for chain in ('LH', 'L', 'H'):
            for numbering in ('Kabat', 'Chothia', 'Martin'):
                endpoint = f"{self.url}/{chain}_Combined_{numbering}/"
                if self.verbose:
                    print(endpoint)
                res = requests.get(endpoint)
                pdb_ids = self.parse_pdb_id(res)
                records = [(pdb_id, chain, numbering) for pdb_id in pdb_ids]
                abYbank_pdb.extend(records)
        abYbank_pdb = list(set(abYbank_pdb))
        return abYbank_pdb

    def parse_pdb_id(self, res):
        soup = BeautifulSoup(res.text, 'html.parser')
        pdb_ids = []
        for line in soup.find_all('a'):
            pdb_id = re.findall(r'>(.*)_\d.pdb', str(line))
            if pdb_id:
                pdb_ids.extend(pdb_id)
        return pdb_ids
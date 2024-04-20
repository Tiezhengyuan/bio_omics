'''
Immune epitope database https://www.iedb.org/
'''
from biosequtils import Dir
import json
import os
import pandas as pd
import numpy as np

from .connector.conn_http import ConnHTTP


class IEDB:
    url = "https://www.iedb.org"
    def __init__(self, local_path:str, version:str=None, overwrite:bool=None):
        self.local_path = local_path
        Dir(self.local_path).init_dir()
        self.version = 'v3' if version is None else version
        self.overwrite = overwrite
        self.meta = {
            'local_path': self.local_path,
            'version': self.version,
        }


    def update_meta(self, info:dict):
        infile = os.path.join(self.local_path, 'meta.json')
        if os.path.isfile(infile):
            with open(infile, 'r') as f:
                self.meta = json.load(f)
        if self.version not in self.meta:
            self.meta[self.version] = {}
        self.meta[self.version][info['name']] = info
        if self.meta:
            with open(infile, 'w') as f:
                json.dump(self.meta, f)

    def pull(self, type:str):
        res = {'name': type,}
        # download
        local_file = self.download_csv(type)
        res['download_file'] = local_file

        # convert to json
        if os.path.isfile(local_file):
            df = pd.read_csv(local_file, compression='zip', header=1, sep=',')
            res['records'] = df.shape[0]
            data = df.to_dict(orient='records')
            json_file = os.path.join(self.local_path, f"{type}_{self.version}.json")
            with open(json_file, 'w') as f:
                json.dump(data, f)
            res['json_file'] = json_file
            self.update_meta(res)
        return res

    def csv_to_dict(self, local_csv:str=None):
        if local_csv is None:
            local_csv = self.meta.get('local_csv', '')
        if os.path.isfile(local_csv):
            df = pd.read_csv(local_csv, compression='zip', header=1, sep=',')
            data = df.to_dict(orient='records')
            return data
        return {}
    
    def download_csv(self, type:str):
        '''
        type: epitope|antigen|tcell|bcell|reference etc
        '''
        _names = {
            'antigen': f"antigen_full_{self.version}.zip",
            'epitope': f"epitope_full_{self.version}.zip",
            'tcell': f"tcell_full_{self.version}.zip",
            'bcell': f"bcell_full_{self.version}_single_file.zip",
            'mhc': f"mhc_ligand_full_single_file.zip",
            'reference': f"reference_full_{self.version}.zip",
            'receptor': f"receptor_full_{self.version}.zip",
            'tcr': f"tcr_full_{self.version}.zip",
            'bcr': f"bcr_full_{self.version}.zip",
            'iedb': f"iedb_3d_full.zip",
        }
        if type in _names:
            conn = ConnHTTP(self.url, self.local_path, self.overwrite)
            file_name = _names[type] 
            end_point = f"/downloader.php?file_name=doc/{file_name}"
            _, local_file = conn.download_file(end_point, file_name)
            return local_file
        return None

    def antigen(self, local_csv:str):
        '''
        unique ID: accession in Antigen IRI
        '''
        df = pd.read_csv(local_csv, compression='zip', header=1, sep=',')
        data_list = df.replace({np.nan:None}).to_dict(orient='records')
        entity_data, n = {}, 1
        for data in data_list:
            key = None
            if data["Antigen IRI"]:
                data['accession'] = os.path.basename(data["Antigen IRI"])
                key = data['accession']
            if data["Organism IRI"] and "Taxon" in data["Organism IRI"]:
                taxon_id = os.path.basename(data["Organism IRI"])
                data['taxon_id'] = taxon_id.split('_')[-1]
            if key is None:
                print(data)
            entity_data[key] = data
            n += 1
            # if n == 100:
            #     break
        return entity_data

    def epitope(self, local_csv:str):
        '''
        unique ID: epitope_id in "IEDB IRI"
        '''
        df = pd.read_csv(local_csv, compression='zip', header=1, low_memory=False)
        data_list = df.replace({np.nan:None}).to_dict(orient='records')
        entity_data, n = {}, 1
        for data in data_list:
            # accession
            # if data["Source Molecule IRI"]:
            #     acc = os.path.basename(data["Source Molecule IRI"])
            #     data['accession'] = acc.split('.')[-1]
            if data["Molecule Parent IRI"]:
                data['accession'] = os.path.basename(data["Molecule Parent IRI"])
            # detect unique key
            key = None
            if data["IEDB IRI"]:
                data['epitope_id'] = os.path.basename(data["IEDB IRI"])
                key = data['epitope_id']
            if key is None:
                print(data)
            entity_data[key] = data
            n += 1
            # if n == 5000:
            #     break
        return entity_data
    

    def reference(self, local_csv:str):
        '''
        unique ID: reference_id in IEDB IRI
        '''
        df = pd.read_csv(local_csv, compression='zip', header=1, sep=',')
        entity_data, n = {}, 1
        data_list = df.replace({np.nan:None}).to_dict(orient='records')
        for data in data_list:
            if data["IEDB IRI"]:
                data['reference_id'] = str(os.path.basename(data["IEDB IRI"]))
                key = data['reference_id']
                entity_data[key] = data
        return entity_data
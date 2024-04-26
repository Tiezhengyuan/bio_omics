'''
Immune epitope database https://www.iedb.org/
'''
from biosequtils import Dir
import json
import os
import pandas as pd
import numpy as np

from ..connector.conn_http import ConnHTTP


class IEDB:
    url = "https://www.iedb.org"
    source = "IEDB"
    meta_file_name = 'IEDB_meta.json'

    def __init__(self, local_path:str, version:str=None, overwrite:bool=None):
        '''
        args: mode: normal(default), debug, update
        '''
        self.local_path = local_path
        Dir(self.local_path).init_dir()
        self.version = 'v3' if version is None else version
        self.overwrite = overwrite
        self.get_meta()
    
    def get_meta(self):
        self.meta_file = os.path.join(self.local_path, self.meta_file_name)
        if os.path.isfile(self.meta_file):
            print('get meta.')
            with open(self.meta_file, 'r') as f:
                self.meta = json.load(f)
        else:
            self.meta = {
                'local_path': self.local_path,
                'source': self.source,
                'version': self.version,
            }
    
    def save_meta(self, local_path:str=None):
        outfile = os.path.join(local_path, self.meta_file_name) \
            if local_path else self.meta_file
        with open(outfile, 'w') as f:
            json.dump(self.meta, f, indent=4)
        return outfile

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

    def csv_to_dict(self, local_csv):
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

    def download_json(self, type:str):
        '''
        type: epitope|antigen|tcell|bcell|reference etc
        '''
        _names = {
            'antigen': f"antigen_full_{self.version}_json.zip",
            'epitope': f"epitope_full_{self.version}_json.zip",
            'tcell': f"tcell_full_{self.version}_json.zip",
            'bcell': f"bcell_full_{self.version}_json.zip",
            'mhc': "mhc_ligand_full_json.zip",
            'reference': f"reference_full_{self.version}_json.zip",
            'tcr': f"tcr_full_{self.version}_json.zip",
            'bcr': f"bcr_full_{self.version}_json.zip",
        }
        if type in _names:
            conn = ConnHTTP(self.url, self.local_path, self.overwrite)
            file_name = _names[type] 
            end_point = f"/downloader.php?file_name=doc/{file_name}"
            _, local_file = conn.download_file(end_point, file_name)
            return local_file
        return None

    def read_json(self, local_file:str) -> list:
        # read json
        data_list = []
        from zipfile import ZipFile
        with ZipFile(local_file) as zf:
            for file in zf.namelist():
                with zf.open(file) as f:
                    data = json.load(f)
                    data_list += data['Data']
        return data_list
    
    def antigen_csv(self) -> dict:
        # download
        file_name = f"antigen_full_{self.version}.zip"
        conn = ConnHTTP(self.url, self.local_path, self.overwrite)
        end_point = f"/downloader.php?file_name=doc/{file_name}"
        _, local_file = conn.download_file(end_point, file_name)
        # convert
        df = pd.read_csv(local_file, compression='zip', header=1, sep=',')
        data_list = df.replace({np.nan:None}).to_dict(orient='records')
        entity_data = self.convert_data(data_list, 'accession')
        return entity_data

    def antigen_json(self) -> dict:
        '''
        unique ID: accession number
        '''
        # download
        file_name = f"antigen_full_{self.version}_json.zip"
        conn = ConnHTTP(self.url, self.local_path, self.overwrite)
        end_point = f"/downloader.php?file_name=doc/{file_name}"
        _, local_file = conn.download_file(end_point, file_name)
       
        # retrieve data
        data_list = self.read_json(local_file)
        entity_data = self.convert_data(data_list, 'accession')
        return entity_data


    def epitope_csv(self) -> dict:
        '''
        unique ID: epitope_id in "IEDB IRI"
        '''
        # download
        file_name = f"epitope_full_{self.version}.zip"
        conn = ConnHTTP(self.url, self.local_path, self.overwrite)
        end_point = f"/downloader.php?file_name=doc/{file_name}"
        _, local_file = conn.download_file(end_point, file_name)
        # convert
        df = pd.read_csv(local_file, compression='zip', header=1, low_memory=False)
        data_list = df.replace({np.nan:None}).to_dict(orient='records')
        entity_data = self.convert_data(data_list, 'epitope_id')
        return entity_data

    def epitope_json(self) -> dict:
        '''
        '''
        # download
        file_name = f"epitope_full_{self.version}_json.zip"
        conn = ConnHTTP(self.url, self.local_path, self.overwrite)
        end_point = f"/downloader.php?file_name=doc/{file_name}"
        _, local_file = conn.download_file(end_point, file_name)
        # retrieve data
        data_list = self.read_json(local_file)
        entity_data = self.convert_data(data_list, 'epitope_id')
        return entity_data
        
    def tcr_csv(self):
        '''
        T-cell receptor
        unique ID: receptor_id
        '''
        # download
        file_name = f"tcr_full_{self.version}.zip"
        conn = ConnHTTP(self.url, self.local_path, self.overwrite)
        end_point = f"/downloader.php?file_name=doc/{file_name}"
        _, local_file = conn.download_file(end_point, file_name)
        # convert
        df = pd.read_csv(local_file, compression='zip', header=1, low_memory=False)
        data_list = df.replace({np.nan:None}).to_dict(orient='records')
        entity_data = self.convert_data(data_list, 'receptor_id')
        return entity_data

    def tcr_json(self):
        '''
        T-cell receptor
        unique ID: receptor_id
        '''
        # download
        file_name = f"tcr_full_{self.version}_json.zip"
        conn = ConnHTTP(self.url, self.local_path, self.overwrite)
        end_point = f"/downloader.php?file_name=doc/{file_name}"
        _, local_file = conn.download_file(end_point, file_name)
        # convert
        data_list = self.read_json(local_file)
        entity_data = self.convert_data(data_list, 'receptor_id')
        return entity_data

    def tcell_csv(self):
        '''
        T-cell
        unique ID: assay_id
        '''
        # download
        file_name = f"tcell_full_{self.version}.zip"
        conn = ConnHTTP(self.url, self.local_path, self.overwrite)
        end_point = f"/downloader.php?file_name=doc/{file_name}"
        _, local_file = conn.download_file(end_point, file_name)
        # convert
        df = pd.read_csv(local_file, compression='zip', header=1, low_memory=False)
        data_list = df.replace({np.nan:None}).to_dict(orient='records')
        entity_data = self.convert_data(data_list, 'assay_id')
        return entity_data

    def tcell_json(self):
        '''
        T-cell
        unique ID: assay_id
        '''
        # download
        file_name = f"tcell_full_{self.version}_json.zip"
        conn = ConnHTTP(self.url, self.local_path, self.overwrite)
        end_point = f"/downloader.php?file_name=doc/{file_name}"
        _, local_file = conn.download_file(end_point, file_name)
        # convert
        data_list = self.read_json(local_file)
        entity_data = self.convert_data(data_list, 'assay_id')
        return entity_data
    
    def bcr_csv(self):
        '''
        B-cell receptor
        unique ID: receptor_id
        '''
        # download
        file_name = f"bcr_full_{self.version}.zip"
        conn = ConnHTTP(self.url, self.local_path, self.overwrite)
        end_point = f"/downloader.php?file_name=doc/{file_name}"
        _, local_file = conn.download_file(end_point, file_name)
        # convert
        df = pd.read_csv(local_file, compression='zip', header=1, low_memory=False)
        data_list = df.replace({np.nan:None}).to_dict(orient='records')
        entity_data = self.convert_data(data_list, 'receptor_id')
        return entity_data

    def bcr_json(self):
        '''
        B-cell receptor
        unique ID: receptor_id
        '''
        # download
        file_name = f"bcr_full_{self.version}_json.zip"
        conn = ConnHTTP(self.url, self.local_path, self.overwrite)
        end_point = f"/downloader.php?file_name=doc/{file_name}"
        _, local_file = conn.download_file(end_point, file_name)
        # convert
        data_list = self.read_json(local_file)
        entity_data = self.convert_data(data_list, 'receptor_id')
        return entity_data

    def bcell_csv(self):
        '''
        B-cell
        unique ID: assay_id
        '''
        # download
        file_name = f"bcell_full_{self.version}_single_file.zip"
        conn = ConnHTTP(self.url, self.local_path, self.overwrite)
        end_point = f"/downloader.php?file_name=doc/{file_name}"
        _, local_file = conn.download_file(end_point, file_name)
        # convert
        df = pd.read_csv(local_file, compression='zip', header=1, low_memory=False)
        data_list = df.replace({np.nan:None}).to_dict(orient='records')
        entity_data = self.convert_data(data_list, 'assay_id')
        return entity_data

    def bcell_json(self):
        '''
        B-cell
        unique ID: assay_id
        '''
        # download
        file_name = f"bcell_full_{self.version}_json.zip"
        conn = ConnHTTP(self.url, self.local_path, self.overwrite)
        end_point = f"/downloader.php?file_name=doc/{file_name}"
        _, local_file = conn.download_file(end_point, file_name)
        # convert
        data_list = self.read_json(local_file)
        entity_data = self.convert_data(data_list, 'assay_id')
        return entity_data

    def mhc_csv(self):
        '''
        major histocompatibility complex
        unique ID: assay_id
        '''
        # download
        file_name = f"mhc_ligand_full_single_file.zip"
        conn = ConnHTTP(self.url, self.local_path, self.overwrite)
        end_point = f"/downloader.php?file_name=doc/{file_name}"
        _, local_file = conn.download_file(end_point, file_name)
        # convert
        df = pd.read_csv(local_file, compression='zip', header=1, low_memory=False)
        data_list = df.replace({np.nan:None}).to_dict(orient='records')
        entity_data = self.convert_data(data_list, 'assay_id')
        return entity_data

    def mhc_json(self):
        '''
        major histocompatibility complex
        unique ID: assay_id
        '''
        # download
        file_name = f"mhc_ligand_full_json.zip"
        conn = ConnHTTP(self.url, self.local_path, self.overwrite)
        end_point = f"/downloader.php?file_name=doc/{file_name}"
        _, local_file = conn.download_file(end_point, file_name)
        # convert
        data_list = self.read_json(local_file)
        entity_data = self.convert_data(data_list, 'assay_id')
        return entity_data
    
    def reference_csv(self):
        '''
        referene
        unique ID: reference_id
        '''
        # download
        file_name = f"reference_full_{self.version}.zip"
        conn = ConnHTTP(self.url, self.local_path, self.overwrite)
        end_point = f"/downloader.php?file_name=doc/{file_name}"
        _, local_file = conn.download_file(end_point, file_name)
        # convert
        df = pd.read_csv(local_file, compression='zip', header=1, low_memory=False)
        data_list = df.replace({np.nan:None}).to_dict(orient='records')
        entity_data = self.convert_data(data_list, 'reference_id')
        return entity_data

    def reference_json(self):
        '''
        reference
        unique ID: reference_id
        '''
        # download
        file_name = f"reference_full_{self.version}_json.zip"
        conn = ConnHTTP(self.url, self.local_path, self.overwrite)
        end_point = f"/downloader.php?file_name=doc/{file_name}"
        _, local_file = conn.download_file(end_point, file_name)
        # convert
        data_list = self.read_json(local_file)
        entity_data = self.convert_data(data_list, 'reference_id')
        return entity_data
            
    def convert_data(self, data_list:list, unique_key:str) -> dict:
        '''
        convert list to dict
        detect duplicates and invalid identifier
        '''
        entity_data = {}
        for data in data_list:
            data['source'] = self.source
            ids = self.retrieve_ids(data)
            data.update(ids)
            # detect unique key
            if unique_key in data:
                key = ids[unique_key]
                entity_data[key] = data
            else:
                print(data)
        return entity_data

    def retrieve_ids(self, data:dict) -> dict:
        '''
        retrieve IDs
        '''
        ids = {}
        for key, value in data.items():
            if value:
                value = str(value)
                base_name = os.path.basename(value)
                if "Antigen IRI" in key or 'Molecule Parent IRI' in key:
                    ids['accession'] = base_name
                elif "IEDB IRI" in key:
                    if "reference" in value:
                        ids['reference_id'] = base_name
                    elif "epitope" in value:
                        ids['epitope_id'] = base_name
                    elif "assay" in value:
                        ids['assay_id'] = base_name
                elif "Epitope IRI" in key and "epitope" in value:
                    ids['epitope_id'] = base_name
                elif "Organism" in key and "NCBITaxon" in value:
                    ids['taxon_id'] = base_name.split('_')[-1]
                elif key == 'Group IRI' and 'receptor' in value:
                    ids['receptor_group_id'] = base_name
                elif 'IEDB Receptor ID' in key:
                    ids['receptor_id'] = value
        return ids




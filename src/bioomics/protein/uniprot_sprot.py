"""
FTP of UniProtKB/Swiss-Prot
"""
import os

from ..integrate_data import IntegrateData
from .uniprot import UniProt


class UniProtSprot(UniProt):
    source = "UniProtKB_SwissProt"
    meta_file_name = 'uniprot_sprot_meta.json'

    def __init__(self, local_path:str, overwrite:bool=None):
        super().__init__(local_path, overwrite)
        self.meta['source'] = self.source
    
    def process_epitopes(self, entity_path:str=None):
        '''
        entity: epitope
        '''
        entity_path = entity_path if entity_path \
            else os.path.join(self.local_path, 'epitope')
        self.meta['entity_path'] = entity_path
        self.integrate = IntegrateData(entity_path)
        self.meta = self.integrate.get_meta(self.meta)
        
        dat_gz = self.download_dat()
        parser = self.parse_dat(dat_gz)
        entity_data = self.parse_epitope(parser)
        self.integrate_epitope(entity_data)

        self.integrate.save_meta(self.meta)
        self.integrate.save_index_meta()
        return True
    
    def download_dat(self, unzip:bool=None):
        '''
        download uniprot_sprot.dat.gz  
        '''
        local_file = self.download_file(
            endpoint = '/pub/databases/uniprot/current_release/knowledgebase/complete',
            file_name = 'uniprot_sprot.dat.gz',
            local_path = self.local_path,
            run_gunzip = True if unzip else False,
        )
        return local_file
  
    def integrate_epitope(self, entity_data:dict):
        '''
        integrate eiptope data into json data
        '''
        m, n = 0, 0
        # check if data exists in json
        for json_data in self.integrate.scan():
            acc = json_data.get('key')
            if  acc in entity_data:
                if 'epitopes' not in json_data:
                    json_data['epitopes'] = {}
                json_data['epitopes'][self.source] = entity_data[acc]['epitopes']
                json_data[self.source] = entity_data[acc]['source']
                self.integrate.save_data(json_data)
                del entity_data[acc]
                m += 1
        # export new data
        for acc, data in entity_data.items():
            input = {
                self.source: data['source'],
                'epitopes': {
                    self.source: data['epitopes']
                },
            }
            self.integrate.add_data(input, acc)
            n += 1
        self.meta['updated_epitopes'] = m
        self.meta['new_epitopes'] = n
        self.meta['epitopes'] = n + m

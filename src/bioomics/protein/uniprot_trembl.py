"""
FTP of UniProtKB/Swiss-Prot
"""
import os

from ..integrate_data import IntegrateData
from .uniprot import UniProt

class UniProtTrembl(UniProt):
    source = "UniProt_TrEMBL"
    meta_file_name = 'uniprot_trembl_meta.json'

    def __init__(self, local_path:str, overwrite:bool=None):
        super().__init__(local_path, overwrite)
        self.meta['source'] = self.source
    
    def process_epitopes(self, entity_path:str=None) -> bool:
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
        count = self.integrate_epitope(self.integrate, entity_data)
        self.meta.update(count)

        self.integrate.save_meta(self.meta)
        self.integrate.save_index_meta()
        return True
    
    def download_dat(self, unzip:bool=None) -> str:
        '''
        download uniprot_trembl.dat.gz  
        '''
        local_file = self.download_file(
            endpoint = '/pub/databases/uniprot/current_release/knowledgebase/complete',
            file_name = 'uniprot_trembl.dat.gz',
            local_path = self.local_path,
            run_gunzip = True if unzip else False,
        )
        return local_file
    

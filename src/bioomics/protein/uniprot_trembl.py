"""
FTP of UniProtKB/Swiss-Prot
"""
import os

from biosequtils import Dir
from ..integrate_data import IntegrateData
from .uniprot import UniProt
from .protein_meta import ProteinMeta

class UniProtTrembl(UniProt):
    source = "UniProt_TrEMBL"

    def __init__(self, local_dir:str, result_dir:str, overwrite:bool=None):
        super().__init__(os.path.join(local_dir, self.source), overwrite, False)
        # path: retrieve data
        self.result_dir = result_dir
        self.meta = {
            'url': self.url,
            'local_path': self.local_path,
        }        
    
    def retrieve_epitopes(self, entity_path:str=None) -> bool:
        '''
        entity: epitope
        '''
        # initialize two objects: meta and integrate
        entity = 'epitope'
        pro_meta = ProteinMeta(self.result_dir, entity, self.source)
        self.meta = pro_meta.get_meta(self.meta)
        self.integrate = IntegrateData(self.meta)

        # download data to the local
        dat_gz = self.download_dat()

        # detect proteins annotated with epitopes in uniprot_sprot.dat
        parser = self.parse_dat(dat_gz)
        entity_data = self.parse_epitope(parser)
        self.meta['count'] = self.integrate_epitope(self.integrate, entity_data)

        # save meta
        self.integrate.save_index_meta()
        pro_meta.save_meta(self.meta)
        return True
    
    def download_dat(self) -> str:
        '''
        download uniprot_trembl.dat.gz  
        '''
        local_file = self.download_file(
            endpoint = '/pub/databases/uniprot/current_release/knowledgebase/complete',
            file_name = 'uniprot_trembl.dat.gz',
            local_path = self.local_path,
        )
        return local_file
    

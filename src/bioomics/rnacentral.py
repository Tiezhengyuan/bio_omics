'''
https://rnacentral.org
'''
import os
from .connector.conn_ftp import ConnFTP


class RNACentral(ConnFTP):
    url = "ftp.ebi.ac.uk"
    endpoint = '/pub/databases/RNAcentral/current_release/'

    def __init__(self, local_dir:str, overwrite:bool=None):
        super().__init__(url=self.url, overwrite=overwrite)
        self.local_dir = os.path.join(local_dir, "RNACentral")
    
    def download_sequence(self, file_name:str) -> tuple:
        _endpoint = f'{self.endpoint}sequences/by-database'
        local_file = self.download_file(
            endpoint = _endpoint,
            file_name = file_name,
            local_path = self.local_dir,
            run_gunzip = False,
        )
        return self.local_dir, local_file
'''
https://rnacentral.org
'''
import os

from ..connector.conn_ftp import ConnFTP


class RNACentral(ConnFTP):
    url = "ftp.ebi.ac.uk"
    endpoint = '/pub/databases/RNAcentral/current_release/'

    def __init__(self, local_dir:str, overwrite:bool=None):
        # run_gunzip: the sequence files on the FTP server are gzip compressed
        super().__init__(url=self.url, overwrite=overwrite, run_gunzip=True)
        self.local_dir = os.path.join(local_dir, "RNACentral")
    
    def download_sequence(self, file_name:str) -> tuple:
        # the files are stored as <file_name>.gz and are gunzipped on arrival
        if not file_name.endswith('.gz'):
            file_name = f'{file_name}.gz'
        local_file = self.download_file(
            endpoint = f'{self.endpoint}sequences/by-database',
            file_name = file_name,
            local_path = self.local_dir,
        )
        return self.local_dir, local_file
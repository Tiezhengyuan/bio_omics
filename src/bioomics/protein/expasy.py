"""
FTP of ExPASy: ftp.expasy.org
"""
import os
from ..connector.conn_ftp import ConnFTP

class Expasy(ConnFTP):
    url = "ftp.expasy.org"

    def __init__(self, local_dir:str, overwrite:bool=None):
        super().__init__(url=self.url, overwrite=overwrite)
        self.local_dir = os.path.join(local_dir, 'ExPASy')
    
    def download_data(self):
        '''
        download data including subdirectories and files
        '''
        local_files = self.download_tree(
            local_path = self.local_dir,
            endpoint = 'databases',
        )
        return local_files

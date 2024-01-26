"""
FTP of ExPASy: ftp.expasy.org
"""
import os
from .conn_ftp import ConnFTP

class ConnectExpasy(ConnFTP):

    def __init__(self):
        endpoint = 'ftp.expasy.org'
        super(ConnectExpasy, self).__init__(endpoint)

    
    def download_data(self):
        '''
        download data including subdirectories and files
        '''
        local_files = self.download_tree(
            local_name = 'expasy',
            ftp_path = '/databases'
        )
        return local_files

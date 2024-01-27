"""
download data from NCIB FTP
"""
import os

from .connector.conn_ftp import ConnFTP
from .connector.conn_ftplib import ConnFTPlib

# lock the scope of groups in local version
ANATOMY_GROUPS = ['archaea', 'bacteria', 'fungi', 'invertebrate', 'plant',
    'protozoa', 'vertebrate_mammalian', 'vertebrate_other', 'viral', ]


class NCBI(ConnFTP):
    url = 'ftp.ncbi.nlm.nih.gov'

    def __init__(self, local_dir:str, overwrite:bool=None):
        super().__init__(url=self.url, overwrite=overwrite)
        self.local_dir = os.path.join(local_dir, "NCBI")

    def download_assembly_summary(self, groups:list=None):
        '''
        download assembly_summary.txt. That is genome metadata
        '''
        if not groups:
            groups = ANATOMY_GROUPS
        groups = [i for i in groups if i in ANATOMY_GROUPS]
        
        res = {}
        for antonomy in groups:
            outdir = os.path.join(self.local_dir, 'assembly_summary', antonomy)
            local_file = self.download_file(
                endpoint = f'genomes/refseq/{antonomy}/',
                file_name = 'assembly_summary.txt',
                local_path = outdir
            )
            res[antonomy] = local_file
        return self.local_dir, res

    def download_genome(self, ftp_path:str, specie:str, version:str=None):
        '''
        download genome including subdirectories and files
        '''
        local_path = os.path.join(self.local_dir, 'genome', specie)
        if version:
            local_path = os.path.join(local_path, version)

        # download sequences and annotations
        local_files = self.download_files(
            endpoint = ftp_path.replace('https://ftp.ncbi.nlm.nih.gov/', ''),
            match = '.gz',
            local_path = local_path,
        )
        return local_path, local_files
    
    def download_gene_data(self):
        '''
        download /gene/DATA including subdirectories and files
        '''
        local_files = self.download_tree(
            local_path = os.path.join(self.local_dir, 'gene', 'DATA'),
            endpoint = 'gene/DATA',
            match = '.gz$'
        )
        return local_files

    def download_pubmed(self):
        '''
        download /PubMed including subdirectories and files
        '''
        res = ConnFTPlib.download_tree(
            ftp_endpoint = self.url,
            ftp_path = '/pubmed',
            match = '.gz',
            local_path = self.local_dir
        )
        return res
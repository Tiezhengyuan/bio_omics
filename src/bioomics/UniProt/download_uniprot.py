"""
FTP of UniProtKB/Swiss-Prot
"""
import os
from pathlib import Path

from ..connector.conn_ftp import ConnFTP

class DownloadUniprot(ConnFTP):
    url = "ftp.uniprot.org"
    source = 'UniProt'

    def __init__(self, outdir, overwrite:bool=False, run_gunzip:bool=False):
        # outdir may be a str: Path is required for the '/' join below
        outdir = Path(outdir)
        super().__init__(self.url, overwrite, run_gunzip)
        self.local_path = outdir / self.source
        self.local_path.mkdir(parents=True, exist_ok=True)

    def uniprot_sprot(self):
        '''
        download uniprot_sprot.dat.gz  
        '''
        local_files = []
        pool = [f'uniprot_sprot.{i}.gz' for i in ('dat', 'fasta', 'xml')]
        for file_name in pool:
            local_file = self.local_path / file_name
            print(local_file)
            if self.overwrite or not local_file.is_file():
                local_file = self.download_file(
                    endpoint = '/pub/databases/uniprot/current_release/knowledgebase/complete',
                    file_name = file_name,
                    local_path = self.local_path,
                )
            local_files.append(local_file)
        return local_files
    
    def uniprot_trembl(self) -> list:
        '''
        download uniprot_trembl.dat.gz  
        '''
        local_files = []
        pool = [f'uniprot_trembl.{i}.gz' for i in ('dat', 'fasta', 'xml')]
        for file_name in pool:
            local_file = self.local_path / file_name
            if self.overwrite or not local_file.is_file():
                local_file = self.download_file(
                    endpoint = '/pub/databases/uniprot/current_release/knowledgebase/complete',
                    file_name = file_name,
                    local_path = self.local_path,
                )
            local_files.append(local_file)
        return local_files

    def uniref50(self) -> Path:
        '''
        download uniref50.fasta.gz  
        '''
        # NB: file on FTP is uniref50.fasta.gz (the pre-check used to test
        # the wrong name 'uniref50.gz' and never detected the local copy)
        local_file = self.local_path / 'uniref50.fasta.gz'
        if self.overwrite or not local_file.is_file():
            local_file = self.download_file(
                endpoint = '/pub/databases/uniprot/current_release/uniref/uniref50',
                file_name = 'uniref50.fasta.gz',
                local_path = self.local_path,
            )
        return Path(local_file)

    def uniref100(self) -> Path:
        '''
        download uniref100.fasta.gz  
        '''
        local_file = self.local_path / 'uniref100.fasta.gz'
        if self.overwrite or not local_file.is_file():
            local_file = self.download_file(
                endpoint = '/pub/databases/uniprot/current_release/uniref/uniref100',
                file_name = 'uniref100.fasta.gz',
                local_path = self.local_path,
            )
        return Path(local_file)

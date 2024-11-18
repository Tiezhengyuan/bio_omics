"""
FTP of UniProtKB/Swiss-Prot
"""
from Bio import SeqIO
from biosequtils import Dir
import gzip
import os
from typing import Iterable

from ..connector.conn_ftp import ConnFTP

class UniProt(ConnFTP):
    url = "ftp.uniprot.org"
    source = 'UniProt'

    def __init__(self, local_dir:str, overwrite:bool=None, run_gunzip:bool=None):
        super().__init__(self.url, overwrite, run_gunzip)
        self.local_path = os.path.join(local_dir, self.source)
        Dir(self.local_path).init_dir()
        self.meta = {
            'url': self.url,
            'local_path': self.local_path,
        }   

    def download_uniprot_sprot(self):
        '''
        download uniprot_sprot.dat.gz  
        '''
        local_file = self.download_file(
            endpoint = '/pub/databases/uniprot/current_release/knowledgebase/complete',
            file_name = 'uniprot_sprot.dat.gz',
            local_path = self.local_path,
        )
        return local_file
    
    def download_uniprot_trembl(self) -> str:
        '''
        download uniprot_trembl.dat.gz  
        '''
        local_file = self.download_file(
            endpoint = '/pub/databases/uniprot/current_release/knowledgebase/complete',
            file_name = 'uniprot_trembl.dat.gz',
            local_path = self.local_path,
        )
        return local_file

    def download_uniref50(self) -> str:
        '''
        download uniref50.dat.gz  
        '''
        local_file = self.download_file(
            endpoint = '/pub/databases/uniprot/current_release/uniref/uniref50',
            file_name = 'uniref50.fasta.gz',
            local_path = self.local_path,
        )
        return local_file

    def download_uniref100(self) -> str:
        '''
        download uniref100.dat.gz  
        '''
        local_file = self.download_file(
            endpoint = '/pub/databases/uniprot/current_release/uniref/uniref100',
            file_name = 'uniref100.fasta.gz',
            local_path = self.local_path,
        )
        return local_file

    def parse_dat(self, dat_file:str) -> Iterable:
        '''
        attributes:
            annotations is dict type
            features: attributes are id, qualifiers, location  
            dbxrefs is list type
            id is str
            seq is Sequence object
            other attr: count, description, name, reverse_complement, translate
        '''
        if dat_file.endswith('gz'):
            with gzip.open(dat_file, 'rt') as f:
                for record in SeqIO.parse(f, 'swiss'):
                    yield record
        else:
            with open(dat_file, 'r') as f:
                for record in SeqIO.parse(f, 'swiss'):
                    yield record

    def parse_fasta(self, infile:str, prefix:str) -> Iterable:
        if infile.endswith('gz'):
            with gzip.open(infile, 'rt') as f:
                for record in SeqIO.parse(f, 'fasta'):
                    rec = {
                        "id": record.id.replace(prefix, ''),
                        "name": record.name,
                        "description": record.description,
                        "sequence": str(record.seq),
                    }
                    yield rec
        else:
            with open(infile, 'r') as f:
                for record in SeqIO.parse(f, 'fasta'):
                    rec = {
                        "id": record.id.replace(prefix, ''),
                        "name": record.name,
                        "description": record.description,
                        "sequence": str(record.seq),
                    }
                    yield rec

    def filter_uniref50(self, index_meta:dict):
        '''
        '''
        fa_gz = self.download_uniref50()
        parser = self.parse_fasta(fa_gz, 'UniRef50_')
        
        entity_data, m, n = {}, 0, 0
        for rec in parser:
            acc = rec['id']
            if acc in index_meta:
                entity_data[acc] = {
                    'protein'
                }

        return entity_data
    
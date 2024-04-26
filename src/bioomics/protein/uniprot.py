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

    def __init__(self, local_path:str, overwrite:bool=None):
        super().__init__(url=self.url, overwrite=overwrite)
        self.local_path = self.local_path = os.path.join(local_path, "UniProt")
        Dir(self.local_path).init_dir()
        self.meta = {
            'url': self.url,
            'local_path': local_path,
        }
  
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

    def parse_epitope(self, parser:Iterable):
        '''
        retrieve records according to keywords defined in features
        args: parser is determined by self.parse_dat()
        '''
        data, m, n = {}, 0, 0
        for record in parser:
            for ft in record.features:
                note = ft.qualifiers.get('note', '')
                if 'epitope' in note:
                    if record.id not in data:
                        m += 1
                        refs = dict([tuple(i.split(':', 1)) for i in record.dbxrefs])
                        data[record.id] = {
                            'source': {
                                'accession': record.id,
                                'annotations': record.annotations,
                                'seq': str(record.seq),
                                'dbxrefs': refs,
                            },
                            'epitopes': {},
                        }
                    n += 1
                    data[record.id]['epitopes'][ft.id] = {
                        'id': ft.id,
                        'qualifiers': ft.qualifiers,
                        'seq': str(record.seq[ft.location.start:ft.location.end]),
                    }
        print(m,n)
        return data

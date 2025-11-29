"""
UniProtKB/Swiss-Prot
 - *.dat
 - *.fasta

"""
from Bio import SeqIO
import gzip
import os
from typing import Iterable


class Uniprot:

    def __init__(self, infile:str):
        # *.dat.gz, *.dat ,or fasta
        self.infile = infile

    def dat_iter(self) -> Iterable:
        '''
        attributes:
            annotations is dict type
            features: attributes are id, qualifiers, location  
            dbxrefs is list type
            id is str
            seq is Sequence object
            other attr: count, description, name, reverse_complement, translate
        '''
        dat_file = self.infile
        if dat_file.endswith('gz'):
            with gzip.open(dat_file, 'rt') as f:
                for record in SeqIO.parse(f, 'swiss'):
                    yield record
        else:
            with open(dat_file, 'r') as f:
                for record in SeqIO.parse(f, 'swiss'):
                    yield record

    def filter_dat_iter(self, ids:list):
        record_iter = self.dat_iter()
        for record in record_iter:
            if str(record.id) not in ids:
                yield record

    def xml_iter(self) -> Iterable:
        if self.infile.endswith('gz'):
            with gzip.open(self.infile, 'rt') as f:
                for record in SeqIO.parse(f, 'uniprot-xml'):
                    yield record
        else:
            with open(dat_file, 'r') as f:
                for record in SeqIO.parse(f, 'uniprot-xml'):
                    yield record

    def fasta_iter(self) -> Iterable:
        fa_file = self.infile
        if fa_file.endswith('gz'):
            with gzip.open(fa_file, 'rt') as f:
                for record in SeqIO.parse(f, 'fasta'):
                    yield record
        else:
            with open(fa_file, 'r') as f:
                for record in SeqIO.parse(f, 'fasta'):
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

    def uniref50_iter(self):
        infile = os.path.join(self.local_dir, 'uniref50.fasta.gz')
        with gzip.open(infile, 'rt') as f:
            for record in SeqIO.parse(f, 'fasta'):
                yield rc
    
    def scan_uniref50(self):
        for record in self.uniref50_iter():
            rec = {
                "id": record.id.replace('UniRef50_', ''),
                "name": record.name,
                "description": record.description,
                "sequence": str(record.seq),
            }
            yield rec

    def scan_fa(self, infile:str):
        for record in self.uniref50_iter():
            rec = {
                "id": record.id,
                "name": record.name,
                "description": record.description,
                "sequence": str(record.seq),
                "features": [feat.__dict__ for feat in record.features],
                "letter_annotations": record.letter_annotations,
                "annotations": record.annotations
            }
            yield rec

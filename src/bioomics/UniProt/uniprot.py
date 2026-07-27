"""
UniProtKB/Swiss-Prot
 - *.dat
 - *.fasta

"""
from Bio import SeqIO
import gzip
from pathlib import Path
from typing import Iterable


class Uniprot:

    def __init__(self, uniprot_dir:Path):
        # *.dat.gz, *.dat ,or fasta
        self.uniprot_dir = uniprot_dir

    def dat_iter(self, infile:Path) -> Iterable:
        '''
        attributes:
            annotations is dict type
            features: attributes are id, qualifiers, location  
            dbxrefs is list type
            id is str
            seq is Sequence object
            other attr: count, description, name, reverse_complement, translate
        '''
        if infile.suffix == '.gz':
            with gzip.open(infile, 'rt') as f:
                for record in SeqIO.parse(f, 'swiss'):
                    yield record
        else:
            with open(infile, 'r') as f:
                for record in SeqIO.parse(f, 'swiss'):
                    yield record

    def filter_uniprot_sprot_dat(self, ids:list):
        infile = self.uniprot_dir / 'uniprot_sprot.dat.gz'
        record_iter = self.dat_iter(infile)
        for record in record_iter:
            if str(record.id) not in ids:
                yield record

    def filter_uniprot_trembl_dat(self, ids:list):
        infile = self.uniprot_dir / 'uniprot_trembl.dat.gz'
        record_iter = self.dat_iter(infile)
        for record in record_iter:
            if str(record.id) not in ids:
                yield record

    def xml_iter(self, infile:Path=None) -> Iterable:
        if not infile:
            infile = self.uniprot_dir / 'uniprot_sprot.xml.gz'  
            
        if infile.suffix == '.gz':
            with gzip.open(infile, 'rt') as f:
                for record in SeqIO.parse(f, 'uniprot-xml'):
                    yield record
        else:
            with open(infile, 'r') as f:
                for record in SeqIO.parse(f, 'uniprot-xml'):
                    yield record

    def fasta_iter(self, infile:Path=None) -> Iterable:
        if not infile:
            infile = self.uniprot_dir / 'uniprot_sprot.fasta.gz'
            
        if infile.suffix == '.gz':
            with gzip.open(infile, 'rt') as f:
                for record in SeqIO.parse(f, 'fasta'):
                    yield record
        else:
            with open(infile, 'r') as f:
                for record in SeqIO.parse(f, 'fasta'):
                    yield record

    def parse_fasta(self, prefix:str, infile:Path=None) -> Iterable:
        if not infile:
            infile = self.uniprot_dir / 'uniprot_sprot.fasta.gz'
            
        if infile.surfix == '.gz':
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

    def uniref50_iter(self, infile:Path=None):
        if not infile:
            infile = self.uniprot_dir / 'uniref50.fasta.gz'
            
        with gzip.open(infile, 'rt') as f:
            for record in SeqIO.parse(f, 'fasta'):
                yield record
    
    def scan_uniref50(self):
        for record in self.uniref50_iter():
            rec = {
                "id": record.id.replace('UniRef50_', ''),
                "name": record.name,
                "description": record.description,
                "sequence": str(record.seq),
            }
            yield rec

    def scan_fa(self):
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

from Bio import GenBank, SeqIO
import gzip
from typing import Iterable

class BioHandler:

    @staticmethod
    def parse_gbk(infile) -> Iterable:
        '''
        NB: coerce to str so that pathlib.Path inputs work with .endswith
        '''
        infile = str(infile)
        if infile.endswith('gz'):
            with gzip.open(infile, 'rt') as f:
                for record in GenBank.parse(f):
                    yield record
        else:
            with open(infile, 'r') as f:
                for record in GenBank.parse(f):
                    yield record

    @staticmethod
    def parse_fasta(infile) -> Iterable:
        '''
        NB: coerce to str so that pathlib.Path inputs work with .endswith
        '''
        infile = str(infile)
        if infile.endswith('gz'):
            with gzip.open(infile, 'rt') as f:
                for record in SeqIO.parse(f, 'fasta'):
                    yield record
        else:
            with open(infile, 'r') as f:
                for record in SeqIO.parse(f, 'fasta'):
                    yield record

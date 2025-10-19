import gzip
import json
import os
import re
from pprint import pprint
from Bio.PDB import PDBIO, PDBParser, MMCIFParser, Select
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from biosequtils import Dir


class ProcessFasta:
    def __init__(self, fa_file:str, verbose:bool=False):
        self.fa_file = fa_file
        self.verbose = verbose

    def to_fasta(self, records):
        try:
            with open(self.fa_file, 'w') as f:
                SeqIO.write(records, f, 'fasta')
                if self.verbose:
                    print(f"Export to {self.fa_file}.")
            return 'succeed'
        except Exception as e:
            if self.verbose:
                print(f"Failure: Can't write to {self.fa_file}, error={e}")
            return 'fail'
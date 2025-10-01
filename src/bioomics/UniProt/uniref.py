import os
from Bio import SeqIO
import gzip

class Uniref:
    def __init__(self, local_dir:str):
        self.local_dir=local_dir
    
    def scan_uniref50(self):
        infile = os.path.join(self.local_dir, 'uniref50.fasta.gz')
        with gzip.open(infile, 'rt') as f:
            for record in SeqIO.parse(f, 'fasta'):
                rec = {
                    "id": record.id.replace('UniRef50_', ''),
                    "name": record.name,
                    "description": record.description,
                    "sequence": str(record.seq),
                }
                yield rec


    def scan_fa(self, infile:str):
        with gzip.open(infile, 'rt') as f:
            for record in SeqIO.parse(f, 'fasta'):
                rec = {
                    "id": record.id,
                    "name": record.name,
                    "description": record.description,
                    "sequence": str(record.seq),
                    "features": [feat.__dict__ for feat in record.features],
                    "letter_annotations": record.letter_annotations,
                    "annotations": record.annotations
                }
                print(rec)
                break
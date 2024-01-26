'''
miRBase: https://mirbase.org/
'''
import os
from Bio import SeqIO
from biosequtils import Dir
from io import StringIO
import lxml.html as html
from urllib.request import urlopen


class Mirbase:
    url = 'https://mirbase.org/download/CURRENT/'

    def __init__(self, local_dir:str, overwrite:bool=None):
        self.local_dir = os.path.join(local_dir, "miRBase")
        Dir(self.local_dir).init_dir()
        self.overwrite = True if overwrite else False
    
    def download_hairpin(self):
        local_file = os.path.join(self.local_dir, 'hairpin.fa')
        if self.overwrite or not os.path.isfile(local_file):
            self.parse_fasta(f"{self.url}/hairpin.fa", local_file)
        return self.local_dir, local_file

    def download_mature(self):
        local_file = os.path.join(self.local_dir, 'mature.fa')
        if self.overwrite or not os.path.isfile(local_file):
            self.parse_fasta(f"{self.url}/mature.fa", local_file)
        return self.local_dir, local_file

    def parse_fasta(self, url_path, local_file):
        # parse xml
        parsed = html.parse(urlopen(url_path))
        doc = parsed.getroot()
        for br in doc.xpath("*//br"):
            br.tail = "\n" + br.tail if br.tail else "\n"
        # parse fasta
        fa_io = StringIO(doc.text_content())
        records = SeqIO.parse(fa_io, 'fasta')
        with open(local_file, 'w') as f:
            for rec in records:
                SeqIO.write(rec, f, 'fasta')

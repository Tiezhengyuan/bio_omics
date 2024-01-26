'''
miRBase: https://mirbase.org/
'''
import os
from Bio import SeqIO
from biosequtils import Dir
from io import StringIO
import lxml.html as html
from urllib.request import urlopen



class ConnectMirbase:
    url = 'https://mirbase.org/download/CURRENT/'

    def __init__(self, outdir:str):
        self.dir_local = os.path.join(outdir, "miRBase")
        Dir(self.dir_local).init_dir()
    
    def download_data(self, name:str, overwrite:bool):
        meta = {
            'hairpin': ('hairpin.fa', 'miRNA_hairpin'),
            'mature': ('mature.fa', 'miRNA_mature'),
        }
        file_name, rna_type = meta.get(name, (None, None))
        # download
        if file_name and rna_type:
            local_file = os.path.join(self.dir_local, file_name)
            if overwrite or not os.path.isfile(local_file):
                self.parse_fasta(f"{self.url}/{file_name}", local_file)
        return rna_type, self.dir_local, local_file

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




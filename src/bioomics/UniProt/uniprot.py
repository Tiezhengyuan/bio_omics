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
        # may be a file or a directory: coerce to Path so that the '/'
        # operator and the suffix checks below always work
        self.uniprot_dir = Path(uniprot_dir)

    @staticmethod
    def _open_text(infile):
        '''
        open a plain or gzipped text file
        '''
        infile = Path(infile)
        if infile.suffix == '.gz':
            return gzip.open(infile, 'rt')
        return open(infile, 'r')

    def _default_file(self, file_name:str) -> Path:
        '''
        default input file: uniprot_dir itself if it is a file,
        else file_name inside the uniprot_dir directory
        '''
        if self.uniprot_dir.is_file():
            return self.uniprot_dir
        return self.uniprot_dir / file_name

    def dat_iter(self, infile:Path=None) -> Iterable:
        '''
        attributes:
            annotations is dict type
            features: attributes are id, qualifiers, location  
            dbxrefs is list type
            id is str
            seq is Sequence object
            other attr: count, description, name, reverse_complement, translate
        '''
        if not infile:
            infile = self._default_file('uniprot_sprot.dat.gz')
        with Uniprot._open_text(infile) as f:
            for record in SeqIO.parse(f, 'swiss'):
                yield record

    def filter_uniprot_sprot_dat(self, ids:list):
        '''
        yield records whose accession is NOT in ids
        (skip already processed accessions)
        '''
        ids = set(ids)
        infile = self.uniprot_dir / 'uniprot_sprot.dat.gz'
        record_iter = self.dat_iter(infile)
        for record in record_iter:
            if str(record.id) not in ids:
                yield record

    def filter_uniprot_trembl_dat(self, ids:list):
        '''
        yield records whose accession is NOT in ids
        (skip already processed accessions)
        '''
        ids = set(ids)
        infile = self.uniprot_dir / 'uniprot_trembl.dat.gz'
        record_iter = self.dat_iter(infile)
        for record in record_iter:
            if str(record.id) not in ids:
                yield record

    def xml_iter(self, infile:Path=None) -> Iterable:
        if not infile:
            infile = self._default_file('uniprot_sprot.xml.gz')
        with Uniprot._open_text(infile) as f:
            for record in SeqIO.parse(f, 'uniprot-xml'):
                yield record

    def fasta_iter(self, infile:Path=None) -> Iterable:
        if not infile:
            infile = self._default_file('uniprot_sprot.fasta.gz')
        with Uniprot._open_text(infile) as f:
            for record in SeqIO.parse(f, 'fasta'):
                yield record

    def parse_fasta(self, prefix:str, infile:Path=None) -> Iterable:
        if not infile:
            infile = self._default_file('uniprot_sprot.fasta.gz')
        # NB: attribute is suffix (was misspelled 'surfix': AttributeError)
        with Uniprot._open_text(infile) as f:
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
            infile = self._default_file('uniref50.fasta.gz')
        with Uniprot._open_text(infile) as f:
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

"""
FTP of UniProtKB/Swiss-Prot
"""
import json
import os

from .uniprot import Uniprot
from .parse_dat_record import ParseDatRecord

class ParseUniprot:
    def __init__(self, outdir:str):
        self.outdir = outdir

    def parse_imgt(self, outdir, name:str='sprot'):
        n, m = 0, 0
        data = []
        dat_file = os.path.join(self.outdir, f'uniprot_{name}.dat.gz')
        record_iter = Uniprot(dat_file).dat_iter()
        for record in record_iter:
            rec_parser = ParseDatRecord(record)
            n += 1
            if rec_parser.is_imgt():
                rec = rec_parser.basic()
                data.append(rec)
                m += 1
                if m % 100 == 0:
                    print(m, end=',')
        print(f'uniprot records = {n}. IMGT records={m}')

        outfile = os.path.join(outdir, f'uniprot_{name}_imgt.json')
        with open(outfile, 'w') as f:
            json.dump(data, f, indent=4)
        return outfile

    def parse_drug(self, outdir, name:str='sprot'):
        n, m = 0, 0
        data = []
        dat_file = os.path.join(self.outdir, f'uniprot_{name}.dat.gz')
        print('Path: ', dat_file)
        record_iter = Uniprot(dat_file).dat_iter()
        for record in record_iter:
            rec_parser = ParseDatRecord(record)
            n += 1
            if rec_parser.is_drug():
                rec = rec_parser.basic()
                data.append(rec)
                m += 1
                if m % 100 == 0:
                    print(m, end=',')
        print(f'uniprot records = {n}. drug records={m}')

        outfile = os.path.join(outdir, f'uniprot_{name}_drug.json')
        with open(outfile, 'w') as f:
            json.dump(data, f, indent=4)
        return outfile

    def filter_uniref50(self, index_meta:dict):
        '''
        retrieve UniRef50 sequences for the accessions defined in index_meta
        input: uniref50.fasta.gz inside self.outdir
        '''
        parser = Uniprot(self.outdir).scan_uniref50()

        entity_data = {}
        for rec in parser:
            acc = rec['id']
            if acc in index_meta:
                entity_data[acc] = {
                    'protein': rec['sequence'],
                }

        return entity_data



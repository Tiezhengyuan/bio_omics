"""
FTP of UniProtKB/Swiss-Prot
"""
from Bio import SeqIO
from biosequtils import Dir
import gzip
import os
import json
from typing import Iterable

from ..connector.conn_ftp import ConnFTP
from ..bio_dict import BioDict
from ..integrate_data import IntegrateData

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
        print("Try to detect epitopes...")
        data, m, n = {}, 0, 0
        for record in parser:
            for ft in record.features:
                note = ft.qualifiers.get('note', '')
                if 'epitope' in note:
                    if record.id not in data:
                        data[record.id] = {
                            'accession': record.id,
                            'source': BioDict.swiss_source(record),
                            'epitopes': [],
                        }
                        m += 1
                    n += 1
                    # update epitope to data
                    epitope = BioDict.swiss_feature(record, ft)
                    data[record.id]['epitopes'].append(epitope)
                    print(json.dumps(data[record.id], indent=4))
        print(f"proteins={m},epitopes={n}")
        return data

    def integrate_epitope(self, integrate_obj:IntegrateData, entity_data:dict):
        '''
        integrate eiptope data into json data
        '''
        count = {
            'epitopes': 0,
            'epitope_proteins': len(entity_data),
            'updated_proteins': 0,
            'updated_epitopes': 0,
            'new_epitopes': 0,
            'new_proteins': 0,
        }
        # check if data exists in json
        for json_data in integrate_obj.scan():
            acc = json_data.get('key')
            if  acc in entity_data:
                if 'epitopes' not in json_data:
                    json_data['epitopes'] = {}
                json_data['epitopes'][self.source] = entity_data[acc]['epitopes']
                json_data[self.source] = entity_data[acc]['source']
                # print(json.dumps(json_data, indent=4))
                integrate_obj.save_data(json_data, (self.source, 'epitope'))
                count['updated_epitopes'] += len(entity_data[acc]['epitopes'])
                count['updated_proteins'] += 1
                del entity_data[acc]
        # export new data
        for acc, data in entity_data.items():
            input = {
                self.source: data['source'],
                'epitopes': {
                    self.source: data['epitopes']
                },
            }
            integrate_obj.add_data(input, acc, (self.source, 'epitope'))
            count['new_epitopes'] += len(data['epitopes'])
            count['new_proteins'] += 1
        count['epitopes'] = count['updated_epitopes'] + count['new_epitopes']
        return count

import pandas as pd
import json
import os


from .ncbi import NCBI
from ..bio_handler import BioHandler
from ..bio_dict import BioDict
from ..integrate_data import IntegrateData

class Refseq:
    source = 'NCBI_refseq'
    meta_file_name = 'ncbi_refseq_meta.json'

    def __init__(self, local_path:str):
        self.local_path = local_path

    def process_epitope(self, entity_path):
        entity_path = entity_path if entity_path \
            else os.path.join(self.local_path, 'epitope')
        self.meta = {
            'source': self.source,
            'entity_path': entity_path,
        }
        self.integrate = IntegrateData(entity_path)
        self.meta = self.integrate.get_meta(self.meta)

        # download
        local_files = []
        # NOTE: no epitopes detect in refseq
        # local_files = NCBI(self.local_path).download_refseq_gpff()
        # local_files = NCBI(self.local_path).download_refseq_complete_gpff()

        # detect epitope in feature.note
        for gpff_gz in local_files:
            if os.path.isfile(gpff_gz):
                print(f"Try to detect epitopes in the file {gpff_gz}.")
                entity_data = self.parse_epitope(gpff_gz)
                count = self.integrate_epitope(entity_data)
                print(count)
                for k,v in count.items():
                    if k not in self.meta:
                        self.meta[k] = 0
                    self.meta[k] += v

        self.integrate.save_meta(self.meta)
        self.integrate.save_index_meta()
        return True
    
    def parse_epitope(self, local_gpff:str):
        '''
        detect proteins with epitopes defined
        '''
        data = {}
        for record in BioHandler.parse_gbk(local_gpff):
            for ft in record.features:
                if ft.key == 'Region':
                    feature = BioDict.gbk_feature(record, ft)
                    if 'epitope' in feature.get('qualifiers', {}).get('note', ''):
                        acc = record.locus
                        if acc not in data:
                            data[acc] = {
                                'source': BioDict.gbk_source(record),
                                'epitopes': [],
                            }
                        data[acc]['epitopes'].append(feature)
                        print(feature)
        return data

    def integrate_epitope(self, entity_data:dict):
        '''
        integrate entity_data into json data in entity_path
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
        for json_data in self.integrate.scan():
            acc = json_data.get('key')
            if  acc in entity_data:
                if 'epitopes' not in json_data:
                    json_data['epitopes'] = {}
                json_data['epitopes'][self.source] = entity_data[acc]['epitopes']
                json_data[self.source] = entity_data[acc]['source']
                # print(json.dumps(json_data, indent=4))
                self.integrate.save_data(json_data, (self.source, 'epitope'))
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
            self.integrate.add_data(input, acc, (self.source, 'epitope'))
            count['new_epitopes'] += len(data['epitopes'])
            count['new_proteins'] += 1
        count['epitopes'] = count['updated_epitopes'] + count['new_epitopes']
        return count



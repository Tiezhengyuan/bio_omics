
'''
'''
import os
import json
from typing import Iterable

from ..bio_dict import BioDict
from ..integrate_data import IntegrateData
from .uniprot import UniProt
from ..protein_meta import ProteinMeta


class UniprotEpitope(UniProt):
    entity = 'epitope'

    def __init__(self, local_dir:str, result_dir:str, overwrite:bool=None):
        super().__init__(local_dir, overwrite, False)
        self.result_dir = result_dir

    def __call__(self, entity_path:str=None) -> bool:
        '''
        entity: epitope
        '''
        iebd_meta = ProteinMeta(self.result_dir, self.entity, 'IEDB')
        iebd_index_meta = iebd_meta.get_index_meta()
        # update self.meta
        pro_meta = ProteinMeta(self.result_dir, self.entity, self.source)
        self.meta = pro_meta.get_meta(self.meta)

        # initialize objects: integrate
        self.integrate = IntegrateData(self.meta)

        # detect proteins annotated with epitopes
        # uniprot_sprot.dat
        # download data to the local
        print(len(iebd_index_meta))
        dat_gz = self.download_uniprot_sprot()
        parser = self.parse_dat(dat_gz)
        entity_data = self.parse_epitope(parser, iebd_index_meta)
        print(len(iebd_index_meta))
        data_source = (self.source, 'uniprot_sprot', self.entity)
        self.meta['uniprot_sprot'] = self.integrate_epitope(entity_data, data_source)

        dat_gz = self.download_uniprot_trembl()
        parser = self.parse_dat(dat_gz)
        entity_data = self.parse_epitope(parser, iebd_index_meta)
        print(len(iebd_index_meta))
        data_source = (self.source, 'uniprot_trembl', self.entity)
        self.meta['uniprot_trembl'] = self.integrate_epitope(entity_data, data_source)

        # retrieve amino acid for IEDB
        # self.retrieve_protein()

        # save meta
        self.integrate.save_index_meta()
        pro_meta.save_meta(self.meta)
        return True


    def parse_epitope(self, parser:Iterable, iebd_index_meta:dict):
        '''
        retrieve records according to keywords defined in features
        args: parser is determined by self.parse_dat()
        '''
        print("Try to detect epitopes...")
        entity_data, m, n, p = {}, 0, 0, 0
        for record in parser:
            tag = 0
            for ft in record.features:
                note = ft.qualifiers.get('note', '')
                if 'epitope' in note:
                    tag = 1
                    if record.id not in entity_data:
                        entity_data[record.id] = {
                            'accession': record.id,
                            'source': BioDict.swiss_source(record),
                            'epitopes': [],
                        }
                        m += 1
                    n += 1
                    # update epitope to data
                    epitope = BioDict.swiss_feature(record, ft)
                    entity_data[record.id]['epitopes'].append(epitope)
                    # print(json.dumps(entity_data[record.id], indent=4))
            if tag == 0 and record.id in iebd_index_meta:
                entity_data[record.id] = {
                    'accession': record.id,
                    'source': BioDict.swiss_source(record),
                    'epitopes': [],
                }
                p += 1
                del iebd_index_meta[record.id]

        print(f"proteins={m}, epitopes={n}, proteins in IEDB={p}")
        self.meta[self.source] = {
            'proteins_with_epitopes': m,
            'epitopes': n,
            'proteins_in_IEDB': p,
        }
        return entity_data

    def integrate_epitope(self, entity_data:dict, data_source:tuple):
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
        for json_data in self.integrate.scan():
            acc = json_data.get('key')
            if  acc in entity_data:
                # key: 'epitopes'
                if 'epitopes' not in json_data:
                    json_data['epitopes'] = {}
                json_data['epitopes'][self.source] = entity_data[acc]['epitopes']
                # key: UniProt
                json_data[self.source] = entity_data[acc]['source']
                # print(json.dumps(json_data, indent=4))
                self.integrate.save_data(json_data, data_source)
                count['updated_epitopes'] += len(entity_data[acc]['epitopes'])
                count['updated_proteins'] += 1
                del entity_data[acc]
        # export new data
        for acc, data in entity_data.items():
            input = {
                self.source: data['source'],
                'epitopes': {
                    self.source: data.get('epitopes', []),
                },
            }
            self.integrate.add_data(input, acc, data_source)
            count['new_epitopes'] += len(data['epitopes'])
            count['new_proteins'] += 1
        count['epitopes'] = count['updated_epitopes'] + count['new_epitopes']
        return count

    def retrieve_protein(self):
        '''
        retrieve protein aa seq based on accession defined by IEDB
        '''
        # get protein accessions defined in IEDB
        # iedb_index_file = os.path.join(self.meta['meta_path'], 'IEDB_index.json')
        # if os.path.isfile(iedb_index_file):
        #     with open(iedb_index_file, 'r') as f:
        #         index_meta = json.load(f)
        #         pro_acc = list(index_meta)
        # print(pro_acc)

        fa_gz = self.download_uniref50()
        print(fa_gz)
        parser = self.parse_fasta(fa_gz, 'UniRef50_')
        for rec in parser:
            print(rec)
            break
        pass
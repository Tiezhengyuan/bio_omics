'''
build data mapping of epitopes from IEDB
entity is proteins with epitopes
'''
from biosequtils import Dir
import os

from .iedb import IEDB
from ..integrate_data import IntegrateData

class IEDBEpitope(IEDB):
    key = 'accession'

    def __init__(self, local_path:str, entity_path:str=None):
        super().__init__(local_path, None, False)
        self.meta['entity'] = 'epitope'
        self.meta['entity_path'] = entity_path if entity_path \
            else os.path.join(self.meta['local_path'], 'epitope')
        Dir(self.meta['entity_path']).init_dir()
        self.integrate = None
    
    def process(self):
        self.integrate = IntegrateData(self.meta['entity_path'])

        print("Process epitopes")
        self.integrate_epitope()
        print("Integrate antigen")
        self.integrate_antigen()
        print("Integrate MHC")
        entity_data = self.mhc_json()
        self.integrate_epitope_related(entity_data, 'MHC')
        print("Integrate B-cell")
        entity_data = self.bcell_json()
        self.integrate_epitope_related(entity_data, 'b_cell')
        print("Integrate B-cell receptor")
        entity_data = self.bcr_json()
        self.integrate_epitope_related(entity_data, 'b_cell_receptor')
        print("Integrate T-cell")
        entity_data = self.tcell_json()
        self.integrate_epitope_related(entity_data, 't_cell')
        print("Integrate T-cell receptor")
        entity_data = self.tcr_json()
        self.integrate_epitope_related(entity_data, 't_cell_receptor')

        # 
        self.integrate.save_index_meta()
        self.save_meta(self.meta['entity_path'])
        return True

    def integrate_epitope(self):
        '''
        epitopes are organized by protein identified by accession
        '''
        entity_data = self.epitope_json()
        # aggregate epitopes by protein
        agg, m, n = {}, 0, 0
        for epitope_id, data in entity_data.items():
            m += 1
            if self.key in data:
                acc = data[self.key]
                if acc not in agg:
                    agg[acc] = {}
                    n += 1
                agg[acc][epitope_id] = data
        self.meta['epitopes'] = m
        self.meta['proteins'] = n
        # check if data exists in json
        for json_data in self.integrate.scan():
            acc = json_data.get('key')
            if  acc in agg:
                if 'epitopes' not in json_data:
                    json_data['epitopes'] = {}
                json_data['epitopes'][self.source] = agg[acc]
                self.integrate.save_data(json_data)
                del agg[acc]
        # export new data
        for acc, data in agg.items():
            input = {
                'epitopes': {self.source: data},
            }
            self.integrate.add_data(input, acc)

    def integrate_antigen(self):
        '''
        Add antigen into eiptope-proteins parsed by accession
        '''
        n = 0
        entity_data = self.antigen_json()
        for json_data in self.integrate.scan():
            acc = json_data.get('key', '')
            if  acc in entity_data:
                if 'antigen' not in json_data:
                    json_data['antigen'] = {}
                json_data['antigen'][self.source] = entity_data[acc]
                self.integrate.save_data(json_data)
                del entity_data[acc]
                n += 1
        self.meta['antigens'] = n
        self.meta['unparsed_antigens'] = len(entity_data)

    def integrate_epitope_related(self, entity_data:dict, inner_key:str):
        # aggregate related data by epitope_id
        agg, n = {}, 0
        for assay_id, data in entity_data.items():
            n += 1
            if 'epitope_id' in data:
                epitope_id = data['epitope_id']
                if epitope_id not in agg:
                    agg[epitope_id] = {}
                agg[epitope_id][assay_id] = data
            # print(json.dumps(data, indent=4))
            # break
        self.meta[inner_key] = n
        self.meta[f"{inner_key}_epitopes"] = len(agg)

        # inject related data into inner_key within that epitope
        m, n = 0, 0
        for json_data in self.integrate.scan():
            tag = 0
            # the key "epitopes" is defined by integrate_epitope()
            if 'epitopes' in json_data and self.source in json_data['epitopes']:
                json_epitopes = json_data['epitopes'][self.source]
                for epitope_id, epitope_data in json_epitopes.items():
                    if epitope_id in agg:
                        epitope_data[inner_key] = agg[epitope_id]
                        # print(json.dumps(agg[epitope_id], indent=4))
                        m += 1
                        n += len(agg[epitope_id])
                        tag = 1
            if tag == 1:
                self.integrate.save_data(json_data)
        self.meta[f"parsed_{inner_key}_epitopes"] = m
        self.meta[f"parsed_{inner_key}"] = n
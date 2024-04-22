'''
build data mapping of antigen from IEDB
'''
from biosequtils import Dir
import os
import json
from .iedb import IEDB
from .integrate_data import IntegrateData

class IEDBEpitope(IEDB):
    key = 'accession'

    def __init__(self, local_path:str):
        super().__init__(local_path, None, False)
        self.meta['entity'] = 'epitope'
        self.meta['entity_path'] = os.path.join(self.meta['local_path'], 'epitope')
        Dir(self.meta['entity_path']).init_dir()
        self.integrate = None
    
    def process(self):
        self.integrate = IntegrateData(self.meta['entity_path'])

        #epitope
        self.integrate_epitope()
        # antigen
        self.integrate_antigen()

        # 
        self.integrate.save_index_meta()
        self.save_meta(self.meta['entity_path'])
        return True

    def integrate_epitope(self):
        entity_data = self.epitope_json()
        # aggregate epitopes by protein
        agg, r, p = {},0, 0
        for epitope_id, data in entity_data.items():
            r += 1
            if self.key in data:
                acc = data[self.key]
                if acc not in agg:
                    agg[acc] = {}
                    p += 1
                agg[acc][epitope_id] = data
        self.meta['epitopes'] = r
        self.meta['proteins_of_epitopes'] = p
        # check if data exists in json
        for json_data in self.integrate.scan():
            acc = json_data.get('key')
            if  acc in agg:
                if 'epitope' not in json_data:
                    json_data['epitope'] = {}
                json_data['epitope'][self.source] = agg[acc]
                self.integrate.save_data(json_data)
                del agg[acc]
        # export new data
        for acc, data in agg.items():
            input = {
                'epitope': {self.source: data},
            }
            self.integrate.add_data(input, acc)

    def integrate_antigen(self):
        entity_data = self.antigen_json()
        # add antigen into proteins with epitopes
        for json_data in self.integrate.scan():
            acc = json_data.get('key', '')
            if  acc in entity_data:
                if 'antigen' not in json_data:
                    json_data['antigen'] = {}
                json_data['antigen'][self.source] = entity_data[acc]
                self.integrate.save_data(json_data)
                del entity_data[acc]


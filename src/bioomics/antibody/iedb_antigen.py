'''
build data mapping of antigen from IEDB

antigen --- eptope --- BCR/TCR/BCell/TCell/MHC
'''
from biosequtils import Dir
import os
import json
from .iedb import IEDB
from ..integrate_data import IntegrateData

class IEDBAntigen(IEDB):
    key = 'accession'

    def __init__(self, local_path:str):
        super().__init__(local_path, None, False)
        self.meta['entity'] = 'antigen'
        self.meta['entity_path'] = os.path.join(self.meta['local_path'], 'antigen')
        Dir(self.meta['entity_path']).init_dir()
        self.integrate = None
    
    def process(self):
        self.integrate = IntegrateData(self.meta['entity_path'])
        
        # antigen
        self.integrate_antigen()
        # epitope
        self.integrate_epitope()
        # MHC
        entity_data = self.mhc_json()
        self.integrate_epitope_related(entity_data, 'MHC')
        # B-cell
        entity_data = self.bcell_json()
        self.integrate_epitope_related(entity_data, 'b_cell')

        # 
        self.integrate.save_index_meta()
        self.save_meta(self.meta['entity_path'])
        return True

    def integrate_antigen(self):
        entity_data = self.antigen_json()
        # check if data exists in json
        for json_data in self.integrate.scan():
            acc = json_data.get('key')
            if  acc in entity_data:
                if 'antigen' not in json_data:
                    json_data['antigen'] = {}
                json_data['antigen'][self.source]=entity_data[acc]
                self.integrate.save_data(json_data)
                del entity_data[acc]
        # add new
        for data in entity_data.values():
            input = {
                'antigen': {self.source: data},
            }
            self.integrate.add_data(input, data.get(self.key))

    def integrate_epitope(self):
        entity_data = self.epitope_json()
        # aggregate epitopes by protein
        agg, m = {},0
        for epitope_id, data in entity_data.items():
            m += 1
            if self.key in data:
                acc = data[self.key]
                if acc not in agg:
                    agg[acc] = {}
                agg[acc][epitope_id] = data

        n=0
        for json_data in self.integrate.scan():
            acc = json_data['key']
            if acc in agg:
                if 'epitopes' not in json_data:
                    json_data['epitopes'] = {}
                json_data['epitopes'][self.source] = agg[acc]
                self.integrate.save_data(json_data)
                n += len(agg[acc])
        print(f"{m}-{n}")


    def integrate_epitope_related(self, entity_data:dict, inner_key:str):
        # aggregate related data by epitope_id
        agg = {}
        for assay_id, data in entity_data.items():
            if 'epitope_id' in data:
                epitope_id = data['epitope_id']
                if epitope_id not in agg:
                    agg[epitope_id] = {}
                agg[epitope_id][assay_id] = data
        # inject related data into inner_key within that epitope
        for json_data in self.integrate.scan():
            tag = 0
            if 'epitopes' in json_data and self.source in json_data['epitopes']:
                json_epitopes = json_data['epitopes'][self.source]
                for epitope_id, epitope_data in json_epitopes.items():
                    if epitope_id in agg:
                        epitope_data[inner_key] = agg[epitope_id]
                        tag = 1
            if tag == 1:
                self.integrate.save_data(json_data)

    

    def _test(self):
        with open('/home/yuan/Downloads/epitope_full_v3.csv', 'r') as f:
            for line in f:
                if "A3X8Q8" in line:
                    print(line)



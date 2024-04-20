'''
build data mapping of antigen from IEDB
'''
from biosequtils import Dir
import os
import json
from .iedb import IEDB
from .integrate_data import IntegrateData

class IEDBEpitope(IEDB):
    entity = 'epitope'
    key = 'accession'

    def __init__(self, local_path:str):
        super().__init__(local_path, None, False)
        self.meta['entity'] = self.entity
        self.meta['entity_path'] = os.path.join(self.meta['local_path'], self.entity)
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


    def download(self, entity:str, file_type:str):
        if file_type == 'csv':
            local_file = self.download_csv(entity)
            key = f"{entity}_csv"
            self.meta[key] = local_file
            return local_file

    def integrate_epitope(self):
        local_file = self.download('epitope', 'csv')
        entity_data = self.epitope(local_file)
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
            acc = json_data.get(self.key)
            if  acc in agg:
                if self.source not in json_data:
                    json_data[self.source] = {}
                json_data[self.source]['epitopes'] = agg[acc]
                self.integrate.save_data(json_data)
                del agg[acc]
        # export new data
        for acc, data in agg.items():
            input = {
                self.source: {
                    'epitopes': data
                }
            }
            self.integrate.add_data(input, acc)


    def integrate_antigen(self):
        self.download('antigen', 'csv')
        entity_data = self.antigen(self.meta['antigen_csv'])
        # add antigen into proteins with epitopes
        for json_data in self.integrate.scan():
            acc = json_data.get('key', '')
            if  acc in entity_data:
                if self.source not in json_data:
                    json_data[self.source] = {}
                json_data[self.source]['antigen'] = entity_data[acc]
                self.integrate.save_data(json_data)
                del entity_data[acc]

    def _test(self):
        import zipfile
        with open('/home/yuan/Downloads/epitope_full_v3.csv', 'r') as f:
            for line in f:
                if "A3X8Q8" in line:
                    print(line)



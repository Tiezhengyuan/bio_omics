'''
build data mapping of antigen from IEDB
'''
from biosequtils import Dir
import os
import json
from .iedb import IEDB
from .integrate_data import IntegrateData

class IEDBAntigen(IEDB):
    source = 'IEDB'
    entity = 'antigen'
    key = 'accession'

    def __init__(self, local_path:str):
        super().__init__(local_path, None, False)
        self.meta['entity'] = self.entity
        self.meta['entity_path'] = os.path.join(self.meta['local_path'], self.entity)
        Dir(self.meta['entity_path']).init_dir()
        self.integrate = None
    
    def process(self):
        self.integrate = IntegrateData(self.meta['entity_path'])
        
        # antigen
        self.integrate_antigen()
        #epitope
        self.integrate_epitope()

        # 
        self.integrate.save_index_meta()
        return True


    def download(self, entity:str, file_type:str):
        if file_type == 'csv':
            local_file = self.download_csv(entity)
            key = f"{entity}_csv"
            self.meta[key] = local_file

    def integrate_antigen(self):
        self.download(self.entity, 'csv')
        entity_data = self.antigen(self.meta['antigen_csv'])
        # check if data exists in json
        for json_data in self.integrate.scan():
            acc = json_data.get(self.key)
            if  acc in entity_data:
                json_data[self.source]=entity_data[acc]
                self.integrate.save_data(json_data)
                del entity_data[acc]
        # add new
        for data in entity_data.values():
            self.integrate.add_data({self.source: data}, data.get(self.key))

    def integrate_epitope(self):
        self.download('epitope', 'csv')
        entity_data = self.epitope(self.meta['epitope_csv'])
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
                json_data[self.source]['epitopes'] = agg[acc]
                self.integrate.save_data(json_data)
                n += len(agg[acc])
            else:
                print(json.dumps(json_data, indent=4))
        print(f"{m}-{n}")

    def _test(self):
        import zipfile
        with open('/home/yuan/Downloads/epitope_full_v3.csv', 'r') as f:
            for line in f:
                if "A3X8Q8" in line:
                    print(line)


    # def integrate_reference(self):
    #     self.download('reference', 'csv')
    #     entity_data = self.reference(self.meta['reference_csv'])
    #     for data in self.integrate.scan():
    #         acc = str(data["# References"])
    #         if acc and acc in entity_data:
    #             data["# References"] = entity_data[acc]
    #             self.integrate.save_data(data)
    #         else:
    #             print(data)

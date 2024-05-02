import pandas as pd
import os
from ..multi_threads import multi_threads

from .ncbi import NCBI
from .parse_id import ParseID
from ..integrate_data import IntegrateData

class IntegrateID:
    source = 'NCBI'
    meta_file_name = 'ncbi_id_meta.json'

    def __init__(self, local_path:str):
        self.local_path = local_path
    
    def ncbi_protein(self, entity_path:str=None):
        '''
        index key is NCBI accession
        '''
        entity_path = entity_path if entity_path \
            else os.path.join(self.local_path, 'ncbi_protein')
        self.meta = {
            'source': self.source,
            'entity_path': entity_path,
            'refseq_uniprotkb': {
                'proteins': 0,
                'updated_proteins': 0,
                'new_proteins': 0,
            },
        }
        self.integrate = IntegrateData(entity_path)
        self.meta = self.integrate.get_meta(self.meta)
        self.index_meta = self.integrate.get_index_meta()

        # donwload and parsing refseq~uniprotkb
        gz_file = NCBI(self.local_path).download_refseq_uniprotkb()
        data_pool_iter = ParseID().gene_refseq_uniprotkb(gz_file, '#NCBI_protein_accession')
        for data_pool in data_pool_iter:
            multi_threads(data_pool, self.integrate_uniprotkb)

        self.integrate.save_index_meta()
        self.integrate.save_meta(self.meta)
        return True

    def integrate_uniprotkb(self, chunk_data):
        '''
        '''
        count = {
            'proteins': 0,
            'updated_proteins': 0,
            'new_proteins': 0,
        }
        print("Try to integrate uniprotkb accession...")
        for acc, group in chunk_data:
            count['proteins'] += 1
            print(acc, end='\t')

            # lift NCBI_tax_id
            tax_arr = group["NCBI_tax_id"].unique()
            acc_source = {
                'NCBI_protein_accession': acc,
                'parent': acc.split('.', 1)[0],
                "NCBI_tax_id": [int(i) for i in tax_arr],
            }
            if tax_arr.all():
                group = group.drop("NCBI_tax_id", axis=1)
            acc_data = group.to_dict(orient="records")

            # check if data exists in json
            related_source, sub_source = 'UniProtKB', f'{self.source}:refseq_uniprotkb'
            if  acc in self.index_meta:
                json_data = self.integrate.get_data(acc)
                json_data.update(acc_source)
                json_data[related_source] = acc_data
                self.integrate.save_data(json_data, self.sub_source)
                count['updated_proteins'] += 1
            # export new data
            else:
                acc_source[related_source] = acc_data
                self.integrate.add_data(acc_source, acc, sub_source)
                count['new_proteins'] += 1
        for k,v in count.items():
            self.meta[sub_source][k] += v
        return count

    def uniprotkb_protein(self, entity_path:str=None):
        '''
        index key is UniProtKB accession
        '''
        entity_path = entity_path if entity_path \
            else os.path.join(self.local_path, 'uniprotkb_protein')
        self.meta = {
            'source': self.source,
            'entity_path': entity_path,
            'uniprotkb_refseq': {
                'proteins': 0,
                'updated_proteins': 0,
                'new_proteins': 0,
            },
        }
        self.integrate = IntegrateData(entity_path)
        self.meta = self.integrate.get_meta(self.meta)
        self.index_meta = self.integrate.get_index_meta()

        # donwload and parsing refseq~uniprotkb
        gz_file = NCBI(self.local_path).download_refseq_uniprotkb()
        data_pool_iter = ParseID().gene_refseq_uniprotkb(gz_file, "UniProtKB_protein_accession")
        for data_pool in data_pool_iter:
            multi_threads(data_pool, self.integrate_refseq)

        self.integrate.save_index_meta()
        self.integrate.save_meta(self.meta)
        return True

    def integrate_refseq(self, chunk_data):
        '''
        '''
        count = {
            'proteins': 0,
            'updated_proteins': 0,
            'new_proteins': 0,
        }
        print("Try to integrate uniprotkb accession...")
        for acc, group in chunk_data:
            count['proteins'] += 1
            print(acc, end='\t')

            # lift NCBI_tax_id
            tax_arr = group["UniProtKB_tax_id"].unique()
            acc_source = {
                "UniProtKB_protein_accession": acc,
                "UniProtKB_tax_id": [int(i) for i in tax_arr],
            }
            if tax_arr.all():
                group = group.drop("UniProtKB_tax_id", axis=1)
            acc_data = group.to_dict(orient="records")

            # check if data exists in json
            related_source, sub_source = 'NCBI_refseq',  f'{self.source}:refseq_uniprotkb'
            if  acc in self.index_meta:
                json_data = self.integrate.get_data(acc)
                json_data.update(acc_source)
                json_data[related_source] = acc_data
                self.integrate.save_data(json_data, sub_source)
                count['updated_proteins'] += 1
            # export new data
            else:
                acc_source[related_source] = acc_data
                self.integrate.add_data(acc_source, acc, sub_source)
                count['new_proteins'] += 1
        for k,v in count.items():
            self.meta[sub_source][k] += v
        return count

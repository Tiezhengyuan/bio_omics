import pandas as pd
import os
from ..concurrency import multi_threads

from .ncbi import NCBI
from .parse_id import ParseID
from .integrate import Integrate
from ..integrate_data import IntegrateData


class IntegrateNCBIProtein(Integrate):
    entity = 'ncbi_protein'

    def __init__(self, local_path:str, overwrite:bool=None, nworks:int=None):
        super().__init__(local_path, overwrite, nworks)

   
    def __call__(self, entity_path:str=None):
        '''
        index key is NCBI accession
        '''
        entity_path = entity_path if entity_path \
            else os.path.join(self.local_path, self.entity)
        self.integrate = IntegrateData(entity_path)
        self.meta = {
            'source': self.source,
            'entity': self.entity,
            'entity_path': entity_path,
        }
        self.meta = self.integrate.get_meta(self.meta)
        self.index_meta = self.integrate.get_index_meta()

        # donwload and parsing refseq~uniprotkb
        self.meta['refseq_uniprotkb'] = dict(DEFAULT_COUNT)
        gz_file = NCBI(self.local_path).download_gene_map('gene_refseq_uniprotkb_collab.gz')
        data_pool_iter = ParseID().gene_map(gz_file, '#NCBI_protein_accession')
        for data_pool in data_pool_iter:
            multi_threads(data_pool, self.integrate_uniprotkb)

        # donwload and parsing gene accession ~ other NCBI accessions
        # self.meta['gene2accession'] = dict(self.default_count)
        # gz_file = NCBI(self.local_path).download_gene_map('gene2accession.gz')
        # data_pool_iter = ParseID().gene_map(gz_file, 'protein_accession.version')
        # for data_pool in data_pool_iter:
        #     multi_threads(data_pool, self.integrate_gene2accession)

        self.integrate.save_index_meta()
        self.integrate.save_meta(self.meta)
        return True

    def integrate_uniprotkb(self, chunk_data):
        '''
        '''
        count = dict(self.default_count)
        print("Try to integrate uniprotkb accession...")
        for acc, group in chunk_data:
            count['proteins'] += 1

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
                if self.overwrite:
                    json_data = self.integrate.get_data(acc)
                    json_data.update(acc_source)
                    json_data[related_source] = acc_data
                    self.integrate.save_data(json_data, sub_source)
                    count['updated_proteins'] += 1
                else:
                    count['skipped_proteins'] += 1
            # export new data
            else:
                acc_source[related_source] = acc_data
                self.integrate.add_data(acc_source, acc, sub_source)
                count['new_proteins'] += 1
        # for k,v in count.items():
        #     self.meta[sub_source][k] += v
        return count

    def integrate_gene2accession(self, chunk_data):
        '''
        '''
        count = dict(DEFAULT_COUNT)
        print("Try to integrate gene2accession...")
        for acc, group in chunk_data:
            count['proteins'] += 1
            acc_data = group.to_dict(orient="records")
            # check if data exists in json
            related_source, sub_source = self.source, f'{self.source}:gene2accession'
            if  acc in self.index_meta:
                if self.overwrite:
                    json_data = self.integrate.get_data(acc)
                    json_data[related_source] = acc_data
                    self.integrate.save_data(json_data, sub_source)
                    count['updated_proteins'] += 1
                else:
                    count['skipped_proteins'] += 1
            # export new data
            else:
                input = {related_source: acc_data}
                self.integrate.add_data(input, acc, sub_source)
                count['new_proteins'] += 1
        for k,v in count.items():
            self.meta[k] += v
        return count
    
    def demo(self, func):
        '''
        for debugging
        '''
        entity_path = os.path.join(self.local_path, 'demo')
        self.integrate = IntegrateData(entity_path)
        self.meta = {
            'source': self.source,
            'entity': self.entity,
            'entity_path': entity_path,
        }
        self.meta = self.integrate.get_meta(self.meta)
        self.index_meta = self.integrate.get_index_meta()

        # donwload and parsing
        self.meta['refseq_uniprotkb'] = dict(DEFAULT_COUNT)
        gz_file = NCBI(self.local_path).download_gene_map('gene_refseq_uniprotkb_collab.gz')
        func(gz_file, self)

        self.integrate.save_index_meta()
        self.integrate.save_meta(self.meta)
        return True
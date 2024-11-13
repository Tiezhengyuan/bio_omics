import pandas as pd
import os

from ..concurrency import *
from .ncbi import NCBI
from .parse_id import ParseID
from .integrate import Integrate
from ..integrate_data import IntegrateData

class IntegrateUniProtKBProtein(Integrate):
    entity = 'uniprotkb_protein'

    def __init__(self, local_path:str, overwrite:bool=None, nworks:int=None):
        super().__init__(local_path, overwrite, nworks)
    
    def __call__(self, entity_path:str=None):
        '''
        index key is UniProtKB accession
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

        parser = ParseID(nworks=self.nworks)
        # donwload and parsing refseq~uniprotkb
        gz_file = NCBI(self.local_path).download_gene_map('gene_refseq_uniprotkb_collab.gz')
        data_pool_iter = parser.gene_map(gz_file, "UniProtKB_protein_accession")
        for data_pool in data_pool_iter:
            results = mp_pool(data_pool, self.integrate_refseq)
            self._join_meta(results, 'ncbi_refseq_uniprotkb')

        self.integrate.save_index_meta()
        self.integrate.save_meta(self.meta)
        return True

    def integrate_refseq(self, chunk_data):
        '''
        '''
        count = dict(self.default_count)
        print("Try to integrate uniprotkb accession...")
        for acc, group in chunk_data:
            count['proteins'] += 1
            existing = self.index_meta.get(acc)
            if existing and self.overwrite is False:
                count['skipped_proteins'] += 1
            else:
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
                if existing:
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
        return count, self.integrate.index_meta

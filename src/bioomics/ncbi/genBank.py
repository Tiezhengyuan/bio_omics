

import json
import os


from .ncbi import NCBI
from ..bio_handler import BioHandler
from ..bio_dict import BioDict
from ..integrate_data import IntegrateData

class GenBank:
    source = 'NCBI_GenBank'
    meta_file_name = 'ncbi_genbank_meta.json'

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
        local_files = [
            '/home/yuan/bio/bio_omics/tests/data/NCBI/protein_fasta/gbbct1.fsa_aa.gz',
            '/home/yuan/bio/bio_omics/tests/data/NCBI/protein_fasta/gbbct2.fsa_aa.gz',
            '/home/yuan/bio/bio_omics/tests/data/NCBI/protein_fasta/gbbct3.fsa_aa.gz', 
        ]
        # NOTE: no epitopes detect in refseq
        # local_files = NCBI(self.local_path).download_protein_fasta()

        # detect epitope in feature.note
        for gz in local_files:
            if os.path.isfile(gz):
                print(f"Try to detect epitopes in the file {gz}.")
                entity_data = self.parse_epitope_proteins(gz)
                # count = self.integrate_epitope(entity_data)
                # print(count)
                # for k,v in count.items():
                #     if k not in self.meta:
                #         self.meta[k] = 0
                #     self.meta[k] += v

        self.integrate.save_meta(self.meta)
        self.integrate.save_index_meta()
        return True

    def parse_epitope_proteins(self, local_gz:str):
        '''
        detect proteins with epitopes defined
        '''
        data = {}
        for record in BioHandler.parse_fasta(local_gz):
            rec = BioDict.fasta(record)
            if 'epitope' in rec['description']:
                acc = rec['accession']
                data[acc] = {
                    'source': rec,
                    'epitopes': [],
                }
                print(json.dumps(data[acc], indent=4))
        return data
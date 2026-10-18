'''
records reday for table insertion
'''
import os
from ..common import ProcessPickle
from ..PDB import Pdb
from ..protein import ParseBlastXml, ParseIgblast
from ..UniProt import ParseDatRecord
from .query_complex import QueryComplex

class BuildRecords:

    def __init__(self, verbose):
        self.verbose = verbose

    
    def uniprot_alphafold_plddt(self, pLDDT_data):
        uniprot_acc = pLDDT_data['uniprot_acc']
        pfile = pLDDT_data['pfile']
        records = []
        for item in pLDDT_data['pLDDT']:
            model_no = item['model'].id
            chain_no = item['chain'].id
            avg_plddt = item['stats'].get('mean_plddt')
            rec = (uniprot_acc, model_no, chain_no, avg_plddt, pfile)
            records.append(rec)
        return records


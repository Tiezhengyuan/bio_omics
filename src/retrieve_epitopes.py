'''
retrieve epitopes from IEDB, UniProt
'''
from bioomics.protein.iedb_epitope import IEDBEpitope
from bioomics.protein.uniprot_sprot import UniProtSprot
from bioomics.protein.uniprot_trembl import UniProtTrembl

def retrieve_epitopes(local_path:str, result_dir:str=None):
    meta = {}

    source = 'IEDB'
    IEDBEpitope(local_path, result_dir).process()
    meta[source] = True

    source = 'UniProtKB_SwissProt'
    UniProtSprot(local_path, result_dir, False).retrieve_epitopes()
    meta[source] = True

    # no eiptope annotation
    # source = 'UniProt_TrEMBL'
    # UniProtTrembl(local_path, result_dir, False).retrieve_epitopes()
    # meta[source] = True

    return meta

local_path = '/home/yuan/data'
result_dir = '/home/yuan/data/omics_data'
retrieve_epitopes(local_path, result_dir)




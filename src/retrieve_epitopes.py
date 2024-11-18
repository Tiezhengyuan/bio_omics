'''
retrieve epitopes from IEDB, UniProt
'''
from bioomics.protein.iedb_epitope import IEDBEpitope
from bioomics.protein.uniprot_epitope import UniprotEpitope

def retrieve_epitopes(local_path:str, result_dir:str=None):
    meta = {}

    # source = 'IEDB'
    # IEDBEpitope(local_path, result_dir).process()
    # meta[source] = True

    source = 'UniProt'
    UniprotEpitope(local_path, result_dir)()
    meta[source] = True

    return meta

local_path = '/home/yuan/data'
result_dir = '/home/yuan/data/omics_data'
retrieve_epitopes(local_path, result_dir)




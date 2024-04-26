'''
retrieve epitopes from IEDB, UniProt
'''
from ..bioomics.immune.iedb_epitope import IEDBEpitope
from ..bioomics.protein.uniprot_sprot import UniProtSprot
from ..bioomics.protein.uniprot_trembl import UniProtTrembl

def retrieve_epitopes(local_path:str, entity_path:str=None):
    meta = {}

    source = 'IEDB'
    try:
        IEDBEpitope(local_path, entity_path).process()
        meta[source] = True
    except Exception as e:
        print(f"Error: Failed in retrieving epitopes from {source}. error={e}")
        meta[source] = False

    source = 'UniProtKB_SwissProt'
    try:
        UniProtSprot(local_path, False).process_epitopes(entity_path)
        meta[source] = True
    except Exception as e:
        print(f"Error: Failed in retrieving epitopes from {source}. error={e}")
        meta[source] = False

    source = 'UniProt_TrEMBL'
    try:
        UniProtTrembl(local_path, False).process_epitopes(entity_path)
        meta[source] = True
    except Exception as e:
        print(f"Error: Failed in retrieving epitopes from {source}. error={e}")
        meta[source] = False

    return meta




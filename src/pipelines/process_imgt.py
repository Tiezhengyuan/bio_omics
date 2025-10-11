'''
example: python src/pipelines/process_imgt.py
functions: sequence alignment
'''
from Bio import SeqIO
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import ProcessBlast, ImgtGenedb, BuildComplex


if __name__ == "__main__":
    args = {
        'data_dir': '/home/yuan/data/pdb/data',
        'collect_dir': '/home/yuan/data/pdb/collect',
        'db_uniprot': '/home/yuan/data/UniProt/db/uniprot_sprot',
        'imgt_dir': '/home/yuan/data/IMGT',
    }
    ig = ImgtGenedb(args['imgt_dir'], True)

    # from gene list to json by spcies
    ig.build_specie_genelist()
    # insert genelist into table imgt_genelist
    genelist = ig.parse_genelist(True)
    records = genelist.to_dict(orient='records')
    BuildComplex(records, True).insert_imgt_genelist()

    # # build fa and db by region
    # ig.build_region_fasta()
    # ig.build_igblastp_region()

    # # build fa and db by organism+spcie
    # ig.build_organism_fasta()
    # ig.build_igblastp_organism()
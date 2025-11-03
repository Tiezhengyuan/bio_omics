'''
example: python src/pipelines/antibody_collect.py
functions:
    - build table imgt_genelist
    - retrieve v-region sequence in fasta
'''
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import ImgtGenedb

def build_genedb_json(args):
    ig = ImgtGenedb(args['imgt_dir'], args['verbose'])
    print('\nRetrieve from gene list to json by spcies')
    ig.build_specie_genelist_json()

def build_region_fasta(args):
    ig = ImgtGenedb(args['imgt_dir'], args['verbose'])

    print("\nRetrieve AA sequences to fasta by v-regions")
    ig.build_region_fasta()
    print("Build igblastp database given region fasta ")
    ig.build_igblastp_region(args['igblast_dir'])

def build_organism_fasta(args):
    ig = ImgtGenedb(args['imgt_dir'], args['verbose'])
    print("\nBuild fa and db by organism+spcie")
    ig.build_organism_fasta()
    print("Build igblastp database given organism fasta ")
    ig.build_igblastp_organism(args['igblast_dir'])


if __name__ == "__main__":
    args = {
        'verbose':True,
        'imgt_dir': '/home/yuan/data/IMGT',
        'igblast_dir': '/home/yuan/data/ncbi-igblast-1.22.0/bin',
    }
    build_genedb_json(args)
    build_region_fasta(args)
    build_organism_fasta(args)

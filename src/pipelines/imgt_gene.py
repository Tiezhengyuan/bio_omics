'''
requirements:
    - IMGT GeneDB
example: python src/pipelines/imgt_gene.py
functions:
    - build table imgt_gene
'''
from collections import defaultdict
from dotenv import load_dotenv
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import ImgtGenedb, BuildComplex

def build_records(args, meta):
    ig = ImgtGenedb(args['imgt_dir'], args['verbose'])
    # read fasta and integrate sequences to self.data
    ig.parse_fasta()

    for specie, sub1 in ig.data.items():
        meta['specie'] += 1
        for gene_name, sub2 in sub1.items():
            for region_name, sub3 in sub2.items():
                for allele_name, sub4 in sub3.items():
                    meta['gene_records'] += 1
                    record = (
                        sub4['record_id'],
                        specie,
                        gene_name[:2].upper(),
                        gene_name[2].upper(),
                        gene_name,
                        allele_name,
                        region_name,
                    )
                    yield record


if __name__ == "__main__":
    meta = defaultdict(int)
    load_dotenv()
    args = {
        'verbose':True,
        'chunk_size': 50,
        'imgt_dir': os.getenv("imgt_dir"),
    }
    record_iter = build_records(args, meta)

    table_name = 'imgt_gene'
    print(f"Try to insert genelist into table {table_name}")
    bc = BuildComplex(record_iter, args['verbose'], args['chunk_size'])
    bc.empty_table(table_name)
    bc.insert_imgt_gene()

    print(meta)
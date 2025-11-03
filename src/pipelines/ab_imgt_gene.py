'''
requirements:
    - IMGT GeneDB
example: python src/pipelines/ab_imgt_gene.py
functions:
    - build table imgt_gene, imgt_geneseq
'''
from collections import defaultdict
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import ImgtGenedb, BuildComplex

def build_records(args):
    ig = ImgtGenedb(args['imgt_dir'], args['verbose'])
    # read fasta and integrate sequences to self.data
    ig.parse_fasta()

    germline_id, records, seq_records = 0, [], []
    for specie, sub1 in ig.data.items():
        for gene_name, sub2 in sub1.items():
            for region_name, sub3 in sub2.items():
                for allele_name, sub4 in sub3.items():
                    germline_id += 1
                    rec = {
                        'germline_id': germline_id,
                        'specie': specie,
                        'isotype': gene_name[:2].upper(),
                        'chain_type': gene_name[2].upper(),
                        'gene_name': gene_name,
                        'allele_name': allele_name,
                        'region_name': region_name,
                    }
                    records.append(rec)
                    for k,v in sub4['seq'].items():
                        seq_rec = {
                            'germline_id': germline_id,
                            'seq_type': k,
                            'seq': v.upper(),
                        }
                        seq_records.append(seq_rec)
    return records, seq_records


if __name__ == "__main__":
    meta = defaultdict(int)
    args = {
        'verbose':True,
        'imgt_dir': '/home/yuan/data/IMGT',
    }
    records, seq_records = build_records(args, meta)
    meta['gene_records'] = len(records)
    meta['geneseq_records'] = len(seq_records)

    table_name = 'imgt_gene'
    print(f"Try to insert genelist into table {table_name}")
    bc = BuildComplex(records, args['verbose'])
    bc.empty_table(table_name)
    bc.insert_imgt_gene()

    table_name = 'imgt_geneseq'
    print(f"Try to insert genelist into table {table_name}")
    bc = BuildComplex(seq_records, args['verbose'])
    bc.empty_table(table_name)
    bc.insert_imgt_geneseq()

    print(meta)
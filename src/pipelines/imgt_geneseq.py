'''
requirements:
    - IMGT GeneDB
example: python src/pipelines/imgt_geneseq.py
functions:
    - build table imgt_geneseq
'''

from helper import *
from bioomics import ImgtGenedb

def build_records(args, meta):
    ig = ImgtGenedb(args['imgt_dir'], args['verbose'])
    # read fasta and integrate sequences to self.data
    ig.parse_fasta()

    for specie, sub1 in ig.data.items():
        meta['specie'] += 1
        for gene_name, sub2 in sub1.items():
            for region_name, sub3 in sub2.items():
                for allele_name, sub4 in sub3.items():
                    for seq_type, seq in sub4['seq'].items():
                        meta['geneseq_records'] += 1
                        record = (
                            sub4['record_id'],
                            seq_type,
                            seq.upper(),
                        )
                        yield record
                        

if __name__ == "__main__":
    args = {
        'verbose':True,
        'chunk_size': 50,
        'imgt_dir': os.getenv("imgt_dir"),
        'table_name': 'imgt_geneseq',
    }
    record_iter = build_records(args, meta)


    print(f"Try to insert genelist into table {args['table_name']}")
    bc = BuildComplex(record_iter, args['verbose'], args['chunk_size'])
    bc.empty_table(args['table_name'])
    bc.insert_batch_records(args['table_name'])

    print(meta)
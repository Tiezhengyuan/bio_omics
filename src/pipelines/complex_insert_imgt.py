'''
example: python src/pipelines/complex_insert_imgt.py
functions:
    - build table imgt_genelist
    - retrieve v-region sequence in fasta
'''
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import ImgtGenedb, BuildComplex

def build_table_imgt_genelist(args):
    ig = ImgtGenedb(args['imgt_dir'], args['verbose'])
    genelist = ig.parse_genelist(to_df=True)
    records = genelist.to_dict(orient='records')

    table_name = 'imgt_genelist'
    print(f"Try to insert genelist into table {table_name}")
    bc = BuildComplex(records, args['verbose'])
    bc.empty_table(table_name)
    bc.insert_imgt_genelist()

def build_table_gene(args):
    ig = ImgtGenedb(args['imgt_dir'], args['verbose'])
    records, seq_records = ig.build_gene_records()

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


if __name__ == "__main__":
    args = {
        'verbose':True,
        'imgt_dir': '/home/yuan/data/IMGT',
    }

    build_table_imgt_genelist(args)
    build_table_gene(args)
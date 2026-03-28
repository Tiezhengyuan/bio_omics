'''
example:
    - parallel -j3 python src/pipelines/clean_files.py 0 ::: {0..6}
    - python src/pipelines/clean_files.py 11
'''
from collections import defaultdict
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import Dir, ProcessFasta, QueryComplex, UpdateComplex


def clean_empty_files(indir):
    n = 0
    for path in Dir(indir).recursive_files():
        if os.path.getsize(path) == 0:
            print(f"Try to delete {path}")
            # os.remove(path)
            n += 1
    print(f"Total {n} files are deleted.")


def clean_table_files(table_name, col_name, file_format, meta):
    print(f"Try to confirm local files defined in {table_name}.{col_name}...")
    files = QueryComplex().unique_values(table_name, col_name)
    n = 0
    for path in files:
        n += 1
        if os.path.isfile(path):
            if file_format == 'fasta':
                res = ProcessFasta(path).clean_fasta()
                if res:
                    meta['having_empty'] += 1
            if os.path.getsize(path) == 0:
                print("The file is empty. Try to delete ", path)
                # delete local file
                os.remove(path)
                # delete row from table
                res = UpdateComplex().delete_rows(table_name, col_name, path)
                meta['deleted'] += res
                meta['no_content'] += 1
            else:
                meta['ready'] += 1
        else:
            print("The file is not existing. Try to delete ", path)
            # delete row from table
            res = UpdateComplex().delete_rows(table_name, col_name, path)
            meta['deleted'] += res
            meta['not_existing'] += 1
        if n % 1000 == 0:
            print(n, end=', ')

if __name__ == "__main__":
    meta = defaultdict(int)
    table_name = sys.argv[1]
    pool = {
        'meta_pdb': ('pdb_file', 'pdb'),
        'meta_chain_faa': ('faa_file', 'fasta'),
        'meta_chain_pdb': ('pdb_file', 'pdb'),
        'meta_combo_pdb': ('pdb_file', 'pdb'),
        'meta_faa': ('faa_file', 'fasta'),
        'meta_pdb_distance': ('pickle_file', 'pickle'),
        'meta_pdb_sasa': ('pickle_file', 'pickle'),
        'uniprot_meta_alphafold': ('pdb_path', 'pdb'),
        'uniprot_alphafold_plddt': ('pickle_file', 'pickle'),
        'meta_blastp5': ('aln_file', 'xml'),
        'meta_pdb_igblastp': ('aln_file', 'text'),
        'meta_bfactor': ('pickle_file', 'pickle'),
        'meta_dabsab': ('summary_file', 'text'),
        'ab_absd_faa': ('faa_file', 'fasta'),
        'ab_absd_aln': ('aln_file', 'text'),
    }
    if table_name in pool:
        col_name, file_format = pool[table_name]
        clean_table_files(table_name, col_name, file_format, meta)
        print('\nSummary:')
        print(f"table:{table_name}, column:{col_name}\n")
        print(meta)
    
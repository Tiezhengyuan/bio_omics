"""
functions:
    - delete empty sub-directories
"""

from src.helper import *


def output_pdb(params, meta):
    for table_name in params['output_pdb_dir'].iterdir():
        path = params['output_pdb_dir'] / table_name
        for root, dirs, files in os.walk(path):
            if dirs == [] and files == []:
                print(root)
                meta[table_name] += 1

if __name__ == "__main__":
    
    output_pdb(params, meta)
    
    footer(meta)
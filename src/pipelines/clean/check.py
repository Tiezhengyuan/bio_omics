
import os

def get_files():
    indir = '/home/yuan/bio/bio_omics/src'
    for root, dirs, files in os.walk(indir):
        for file in files:
            if file.endswith('.py'):
                path = os.path.join(root, file)
                yield path


for path in get_files():
    with open(path, 'r') as f:
        for line in f:
            if 'protein_source' in line:
                print(path)
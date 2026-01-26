'''
example: python src/pipelines/pdb_collect.py
functions: collect data determined by analyze_pdb.py and export to collect/
'''
from Bio import SeqIO
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import Dir, ProcessSeq

def chain_sequences(args):
    proteins, chains = 0, 0
    file_name = 'chain_sequences.faa'
    outfile = os.path.join(args['collect_dir'], file_name)
    with open(outfile, 'w') as OUT:
        file_iter = Dir(args['data_dir']).recursive_files()
        for path in file_iter:
            if path.endswith(file_name):
                proteins += 1
                with open(path, 'r') as IN:
                    parser = SeqIO.parse(IN, 'fasta')
                    for record in parser:
                        SeqIO.write(record, OUT, 'fasta')
                        chains += 1
    print(f"Combined {chains} out of {proteins} files into {outfile}.")

def to_fasta(rows, outfile):
    records = []
    for row in rows:
        _ids = [row['chain_id'], row['uniprot_acc'], row['allele_name']]
        records.append({
            'seq': row['pro_seq'],
            'id': '|'.join(_ids),
        })
    ProcessSeq.to_fasta(records, outfile)

if __name__ == "__main__":
    args = {
        'verbose': True,
        'data_dir': '/home/yuan/data/pdb/data',
        'collect_dir': '/home/yuan/data/pdb/collect',
    }
    # export chain sequences in fasta
    chain_sequences(args)


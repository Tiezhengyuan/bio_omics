'''
example: python src/pipelines/multi_align.py
functions: collect data from tables and perform MSA
'''
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import Dir, QueryComplex, ProcessSeq, AlignSeq

def to_fasta(rows, outfile):
    records = []
    for row in rows:
        _ids = [row['chain_id'], row['uniprot_acc'], row['allele_name']]
        records.append({
            'seq': row['pro_seq'],
            'id': '|'.join(_ids),
        })
    ProcessSeq.to_fasta(records, outfile)

def heavy_chain(args):
    name = 'antibody_heavy'
    # prepare fasta
    q = QueryComplex(args['verbose'])
    rows = q.get_heavy_seq()
    fa_file = os.path.join(args['outdir'], f'{name}.faa')
    to_fasta(rows, fa_file)
    # MSA
    aln_file = os.path.join(args['outdir'], f'{name}.afa')
    AlignSeq(args['verbose']).muscle(fa_file, aln_file, True)

def light_chain(args):
    name = 'antibody_light'
    # prepare fasta
    q = QueryComplex(args['verbose'])
    rows = q.get_light_seq()
    fa_file = os.path.join(args['outdir'], f'{name}.faa')
    to_fasta(rows, fa_file)
    # MSA
    aln_file = os.path.join(args['outdir'], f'{name}.afa')
    AlignSeq(args['verbose']).muscle(fa_file, aln_file, True)

if __name__ == "__main__":
    args = {
        'verbose': True,
        'data_dir': '/home/yuan/data/pdb/data',
        'outdir': '/home/yuan/data/pdb/MSA',
    }
    Dir(args['outdir']).init_dir()
    
    # export chain sequences in fasta and do MSA
    heavy_chain(args)
    light_chain(args)
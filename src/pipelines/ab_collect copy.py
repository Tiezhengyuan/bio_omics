'''
example: python src/pipelines/antibody_collect.py
functions: collect data determined by analyze_pdb.py and export to collect/
'''
import os
import pandas as pd
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import Dir, QueryComplex, ProcessSeq, AlignSeq, ProcessPickle

def to_fasta(rows, outfile):
    records = []
    for row in rows:
        _ids = [row['chain_id'], row['uniprot_acc'], row['allele_name']]
        records.append({
            'seq': row['pro_seq'],
            'id': '|'.join(_ids),
        })
    ProcessSeq.to_fasta(records, outfile)

def collect_heavy_chain(args):
    outdir = args['collect_dir']
    Dir(outdir).init_dir()

    name = 'antibody_heavy'
    # prepare fasta
    q = QueryComplex(args['verbose'])
    rows = q.get_heavy_seq()
    fa_file = os.path.join(outdir, f'{name}.faa')
    to_fasta(rows, fa_file)
    # do MSA
    aln_file = os.path.join(outdir, f'{name}.afa')
    AlignSeq(args['verbose']).muscle(fa_file, aln_file, True)

def collect_light_chain(args):
    outdir = args['collect_dir']
    Dir(outdir).init_dir()

    name = 'antibody_light'
    # prepare fasta
    q = QueryComplex(args['verbose'])
    rows = q.get_light_seq()
    fa_file = os.path.join(outdir, f'{name}.faa')
    to_fasta(rows, fa_file)
    # do MSA
    aln_file = os.path.join(outdir, f'{name}.afa')
    AlignSeq(args['verbose']).muscle(fa_file, aln_file, True)

def collect_msa_distance(args):
    print("\bTry to collect phylotgentic distance based on MSA...")
    res, n = {}, 0
    indir = os.path.join(args['data_dir'], 'data')
    file_iter = Dir(indir).recursive_files()
    for path in file_iter:
        if path.endswith('.habs.pkl'):
            data = ProcessPickle(path).load_pickle()
            chain_type = data['chain_type']
            for method, dist in data['distance'].items():
                key = (chain_type, method)
                if key not in res:
                    res[key] = []
                res[key].append(dist)
            n += 1
    if args['verbose']:
        print(f"Detect {n} *.habs.pk files")

    # export to csv
    collect_dir = os.path.join(args['data_dir'], 'collect')
    for (chain_type, method), dist_pool in res.items():
        df = pd.concat(dist_pool, axis=1).T
        outfile = os.path.join(collect_dir, f'phylo_habs_{chain_type}_{method}.csv')
        df.to_csv(outfile, index_label='chain_id')
        if args['verbose']:
            print(df.shape, outfile)

if __name__ == "__main__":
    args = {
        'verbose': True,
        'data_dir': '/home/yuan/data/pdb',

    }

    # export chain sequences in fasta and do MSA
    # collect_heavy_chain(args)
    # collect_light_chain(args)

    # MSA distance
    collect_msa_distance(args)
'''
example: 
    - parallel -j4 python /home/yuan/bio/bio_omics/src/pipelines/pdb_fasta.py \
        ::: /home/yuan/data/pdb/_list/pdb_*.json
    - python /home/yuan/bio/bio_omics/src/pipelines/pdb_fasta.py \
        /home/yuan/data/pdb/_list/pdb_46.json
functions:
    - retrieve pdb to fasta
    - blastp against uniprot_sprot
'''
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import Pdb, ProcessJson

def to_fasta(args):
    p = Pdb(args['pdb_file'], args['outdir'], False)
    p.load_structure()
    cascade_dir = p.init_cascade_outdir(args['pdb_id'])
    if cascade_dir:
        outfile = os.path.join(cascade_dir, 'chain_sequences.faa')
        if args['overwrite'] or not os.path.isfile(outfile):
            status = p.export_fasta(outfile)
            return status
        else:
            return 'skip'
    return 'null'

def to_chain_fasta(args):
    # export to data/ in default
    p = Pdb(args['pdb_file'], args['outdir'], args['verbose'])
    p.load_structure()
    status = p.export_chains_fasta(args['overwrite'])
    return status

def to_chain_pdb(args):
    # export to data/ in default
    p = Pdb(args['pdb_file'], args['outdir'], args['verbose'])
    p.load_structure()
    status = p.split_by_chain(args['overwrite'])
    return status



if __name__ == "__main__":
    status = ['succeed', 'skip', 'null', 'na', 'fail']
    meta = {
        'to_fasta': {i:0 for i in status},
        'to_chain_fasta': {i:0 for i in status},
        'to_chain_pdb': {i:0 for i in status},
    }
    json_file = sys.argv[1]
    args_pool = ProcessJson.load_json(json_file)
    for args in args_pool:
        args['verbose'] = True
        args['overwrite'] = False
        
        # export chain sequences in fasta
        status = to_fasta(args)
        meta['to_fasta'][status] += 1
        status = to_chain_fasta(args)
        meta['to_chain_fasta'][status] += 1
        status = to_chain_pdb(args)
        meta['to_chain_pdb'][status] += 1
    print(meta)

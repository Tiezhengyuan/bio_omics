'''
example: 
    - parallel -j4 python /home/yuan/bio/bio_omics/src/pipelines/pdb_blast.py \
        ::: /home/yuan/data/pdb/_list/pdb_*.json
    - python /home/yuan/bio/bio_omics/src/pipelines/pdb_blast.py \
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

from bioomics import Pdb, ProcessJson, ProcessBlast


# blastp of chain sequences against uniprot_sprot
def uniprot_blastp(args):
    status = 'null'
    p = Pdb(args['pdb_file'], args['outdir'], args['verbose'])
    cascade_dir = p.init_cascade_outdir(args['pdb_id'])
    query_fa = os.path.join(cascade_dir, 'chain_sequences.faa')
    if os.path.isfile(query_fa):
        for outfmt in (5,):
            aln_file = os.path.join(cascade_dir, f'chain_sequences.uniprot.blastp{outfmt}')
            if args['overwrite'] or not os.path.isfile(aln_file):
                pb = ProcessBlast(verbose=args['verbose'])
                status = pb.blastp(query_fa, args['db_uniprot'], aln_file, outfmt)
            else:
                status = 'skip'
    return status


if __name__ == "__main__":
    status = ['succeed', 'skip', 'null', 'na', 'fail']
    json_file = sys.argv[1]
    meta = {
        'json_file': json_file,
        'blastp_uniprot': {i:0 for i in status},
    }

    args_pool = ProcessJson.load_json(json_file)
    for args in args_pool:
        args['verbose'] = True
        args['overwrite'] = False

        # blastp against uniprot 
        args['db_uniprot'] = '/home/yuan/data/UniProt/db/uniprot_sprot'
        status = uniprot_blastp(args)
        meta['blastp_uniprot'][status] += 1

    print(meta)

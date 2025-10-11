'''
example: 
    - parallel -j4 python src/pipelines/pdb_analyze.py ::: /home/yuan/data/pdb/_list/pdb_*.json
    - python src/pipelines/pdb_analyze.py /home/yuan/data/pdb/_list/pdb_9.json
functions:
    - retrieve pdb to fasta
    - blastp against uniprot_sprot
'''
from Bio import SeqIO
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import PullPdb, Pdb, ProcessJson, ProcessBlast

def to_fasta(args, meta):
    p = Pdb(args['pdb_file'], args['outdir'], False)
    cascade_dir = p.init_cascade_outdir(args['pdb_id'])
    if cascade_dir:
        outfile = os.path.join(cascade_dir, 'chain_sequences.faa')
        if args['overwrite'] or not os.path.isfile(outfile):
            status = p.export_fasta(outfile)
            meta['to_fasta'][status] += 1
        else:
            meta['to_fasta']['skip'] += 1
    else:
        meta['to_fasta']['null'] += 1

def uniprot_blastp(args, meta):
    p = Pdb(args['pdb_file'], args['outdir'], False)
    cascade_dir = p.init_cascade_outdir(args['pdb_id'])
    query_fa = os.path.join(cascade_dir, 'chain_sequences.faa')
    if os.path.isfile(query_fa):
        for outfmt in (5,):
            aln_file = os.path.join(cascade_dir, f'chain_sequences.uniprot.blastp{outfmt}')
            if args['overwrite'] or not os.path.isfile(aln_file):
                pb = ProcessBlast(verbose=False)
                status = pb.blastp(query_fa, args['db_uniprot'], aln_file, outfmt)
                meta['blastp_uniprot'][status] += 1
            else:
                meta['blastp_uniprot']['skip'] += 1
    else:
        meta['blastp_uniprot']['null'] += 1


if __name__ == "__main__":
    json_file = sys.argv[1]
    status = ['succeed', 'skip', 'null', 'na', 'fail']
    meta = {
        'json_file': json_file,
        'to_fasta': {i:0 for i in status},
        'blastp_uniprot': {i:0 for i in status},
    }

    args_pool = ProcessJson.load_json(json_file)
    for args in args_pool:
        # export chain sequences in fasta
        to_fasta(args, meta)
        # blastp against uniprot 
        args['db_uniprot'] = '/home/yuan/data/UniProt/db/uniprot_sprot'
        args['db_vregion'] = '/home/yuan/data/IMGT/GENE-DB/AA-WithoutGaps-F+ORF+inframeP/region_igblastdb/V-REGION'
        uniprot_blastp(args, meta)
    print(meta)

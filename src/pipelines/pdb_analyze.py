'''
example: 
    - parallel -j4 python /home/yuan/bio/bio_omics/src/pipelines/pdb_analyze.py \
        ::: /home/yuan/data/pdb/_list/pdb_*.json
    - python /home/yuan/bio/bio_omics/src/pipelines/pdb_analyze.py \
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

from bioomics import Pdb, ProcessJson, ProcessBlast, ProcessIgblast

def to_fasta(args, meta):
    key = 'to_fasta'
    p = Pdb(args['pdb_file'], args['outdir'], False)
    cascade_dir = p.init_cascade_outdir(args['pdb_id'])
    if cascade_dir:
        outfile = os.path.join(cascade_dir, 'chain_sequences.faa')
        if args['overwrite'] or not os.path.isfile(outfile):
            status = p.export_fasta(outfile)
            if status not in meta['to_fasta']:
                meta[key] = 1
            else:
                meta[key][status] += 1
        else:
            meta[key]['skip'] += 1
    else:
        meta[key]['null'] += 1

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
                if status not in meta['to_fasta']:
                    meta['blastp_uniprot'] = 1
                else:
                    meta['blastp_uniprot'][status] += 1
            else:
                meta['blastp_uniprot']['skip'] += 1
    else:
        meta['blastp_uniprot']['null'] += 1

def vregion_igblastp(args, meta):
    key = 'vregion.igblastp'
    p = Pdb(args['pdb_file'], args['outdir'], args['verbose'])
    cascade_dir = p.init_cascade_outdir(args['pdb_id'])
    query_fa = os.path.join(cascade_dir, 'chain_sequences.faa')
    if os.path.isfile(query_fa):
        outfmt = 3
        aln_file = os.path.join(cascade_dir, f'chain_sequences.{key}{outfmt}')
        if args['overwrite'] or not os.path.isfile(aln_file):
            pb = ProcessIgblast(verbose=args['verbose'])
            status = pb.igblastp(query_fa, args['db_vregion'], aln_file, outfmt)
            if status not in meta['to_fasta']:
                meta[key] = 1
            else:
                meta[key][status] += 1

        else:
            meta[key]['skip'] += 1
    else:
        meta[key]['null'] += 1

def to_chain_fasta(args, meta):
    p = Pdb(args['pdb_file'], args['outdir'], args['verbose'])
    p.export_chains_fasta(args['overwrite'])

def to_chain_pdb(args, meta):
    p = Pdb(args['pdb_file'], args['outdir'], args['verbose'])
    p.split_by_chain(args['overwrite'])


if __name__ == "__main__":
    status = ['succeed', 'skip', 'null', 'na', 'fail']
    json_file = sys.argv[1]
    meta = {
        'json_file': json_file,
        'to_fasta': {i:0 for i in status},
        'blastp_uniprot': {i:0 for i in status},
        'vregion.igblastp': {i:0 for i in status},
    }

    args_pool = ProcessJson.load_json(json_file)
    for args in args_pool:
        args['verbose'] = True
        args['overwrite'] = False
        
        # export chain sequences in fasta
        to_fasta(args, meta)
        to_chain_fasta(args, meta)
        to_chain_pdb(args, meta)
        
        # blastp against uniprot 
        args['db_uniprot'] = '/home/yuan/data/UniProt/db/uniprot_sprot'
        uniprot_blastp(args, meta)
        
        # igblastp agaist IMGT V-region
        args['db_vregion'] = '/home/yuan/data/IMGT/GENE-DB/AA-WithoutGaps-F+ORF+inframeP/region_igblastdb/V-REGION'
        vregion_igblastp(args, meta)
    print(meta)

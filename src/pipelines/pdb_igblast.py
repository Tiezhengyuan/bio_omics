'''
example: 
    - cd /home/yuan/data/ncbi-igblast-1.22.0
    - parallel -j4 python /home/yuan/bio/bio_omics/src/pipelines/pdb_igblast.py \
        ::: /home/yuan/data/pdb/_list/pdb_*.json
    - python /home/yuan/bio/bio_omics/src/pipelines/pdb_igblast.py \
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

from bioomics import Pdb, ProcessJson, ProcessIgblast


def vregion_igblastp(args):
    status = 'null'
    p = Pdb(args['pdb_file'], args['outdir'], args['verbose'])
    cascade_dir = p.init_cascade_outdir(args['pdb_id'])
    query_fa = os.path.join(cascade_dir, 'chain_sequences.faa')
    if os.path.isfile(query_fa):
        outfmt = 3
        aln_file = os.path.join(cascade_dir, f'chain_sequences.vregion.igblastp{outfmt}')
        if args['overwrite'] or not os.path.isfile(aln_file):
            pb = ProcessIgblast(verbose=args['verbose'])
            status = pb.igblastp(query_fa, args['db_vregion'], aln_file, outfmt)
        else:
            status = 'skip'
    return status

if __name__ == "__main__":
    status = ['succeed', 'skip', 'null', 'na', 'fail']
    json_file = sys.argv[1]
    meta = {
        'json_file': json_file,
        'vregion': {i:0 for i in status},
    }

    args_pool = ProcessJson.load_json(json_file)
    for args in args_pool:
        args['verbose'] = True
        args['overwrite'] = False

        # igblastp agaist IMGT V-region
        args['db_vregion'] = '/home/yuan/data/IMGT/GENE-DB/AA-WithoutGaps-F+ORF+inframeP/region_igblastdb/V-REGION'
        status = vregion_igblastp(args)
        meta['vregion'][status] += 1
    print(meta)

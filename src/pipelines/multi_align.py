'''
example: 
    - parallel -j4 python /home/yuan/bio/bio_omics/src/pipelines/multi_align.py \
        ::: /home/yuan/data/pdb/_list/pdb_*.json
    - python /home/yuan/bio/bio_omics/src/pipelines/multi_align.py \
        /home/yuan/data/pdb/_list/pdb_46.json

functions: collect data from tables and perform MSA
'''
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import Dir, QueryComplex, AlignSeq, CalPhylo, \
    ProcessJson, ProcessSeq, ProcessPickle


def align_germline(args:dict, pdb_ids:dict):
    print('Retrieve antibody + germline sequence...')
    qc = QueryComplex(args['verbose'])
    records_pool = qc.germline_human_vregion(pdb_ids.keys())

    print('Do MSA...')
    for chain_type, pool in records_pool.items():
        for item in pool:
            ab, germlines = item
            chain_id = ab['id']
            pdb_id = ab['pdb_id']
            pdb_file = pdb_ids[pdb_id]
            outdir = os.path.join(args['outdir'], pdb_id[:2], pdb_id, 'MSA')
            Dir(outdir).init_dir()
            outprefix = os.path.join(outdir, f"{chain_id}.{chain_type}")

            # build fasta
            fa_file = outprefix + '.habs.faa'
            if not os.path.isfile(fa_file):
                del ab['pdb_id']
                records = [ab,] + germlines
                ProcessSeq.to_fasta(records, fa_file)
            
            # multiple sequence alignment
            aln_file = outprefix + '.habs.afa'
            if not os.path.isfile(aln_file):
                AlignSeq(args['verbose']).muscle(fa_file, aln_file, large_msa=True)

            print('Try to calculate phylogentic distance...')
            if os.path.isfile(aln_file):
                data = CalPhylo().msa_distance(aln_file)
                if data:
                    data.update({
                        'chain_type': chain_type,
                        'chain_id': chain_id,
                        'pdb_id': pdb_id,
                        'fa_file': fa_file,
                        'pdb_file': pdb_file,
                    })
                    pfile = outprefix + '.habs.pkl'
                    ProcessPickle(pfile).save_pickle(data)
                    if args['verbose']:
                        print('Export phylogentic distiance to ', pfile)



##############################
if __name__ == "__main__":
    args = {
        'verbose':True,
        'outdir': '/home/yuan/data/pdb/data',
    }
    json_file = sys.argv[1]
    args_pool = ProcessJson.load_json(json_file)
    pdb_ids = {i['pdb_id']:i['pdb_file'] for i in args_pool}

    align_germline(args, pdb_ids)

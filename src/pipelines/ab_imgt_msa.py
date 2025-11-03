'''
requirements:
    - table pdb_chain, meta_chain_faa
example: 
    - parallel -j4 python src/pipelines/ab_imgt_msa.py ::: {0..46}
    - python src/pipelines/ab_imgt_msa.py 46
functions:
    - collect data from tables and perform MSA
'''
from collections import defaultdict
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import Dir, QueryComplex, AlignSeq, CalPhylo, \
    QueryComplex, ProcessSeq, ProcessPickle

def retrieve_data(pdb_group, args, meta):
    query = f"""
        SELECT G.allele_name AS id, S.seq AS seq
        FROM imgt_gene G, imgt_geneseq S
        WHERE G.chain_type IN {chain_group}
            AND G.germline_id = S.germline_id
            AND G.specie = 'Homo_sapiens'
            AND G.isotype = 'IG'
            AND G.region_name = 'V-REGION'
            AND S.seq_type = 'AA-WithoutGaps-F+ORF+inframeP';
    """
    rows = QueryComplex(args['verbose']).list_data(query)
    return rows


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

            pfile = outprefix + '.habs.pkl'
            if os.path.isfile(aln_file) and (not os.path.isfile(pfile)):
                print('Try to calculate phylogentic distance...')
                data = CalPhylo().msa_distance(aln_file)
                if data:
                    data.update({
                        'chain_type': chain_type,
                        'chain_id': chain_id,
                        'pdb_id': pdb_id,
                        'fa_file': fa_file,
                        'pdb_file': pdb_file,
                    })
                    ProcessPickle(pfile).save_pickle(data)
                    if args['verbose']:
                        print('Export phylogentic distiance to ', pfile)



##############################
if __name__ == "__main__":
    meta = defaultdict(int)
    args = {
        'verbose':True,
        'outdir': '/home/yuan/data/pdb/data',
    }
    pdb_group = sys.argv[1]
    record_iter = retrieve_data(pdb_group, args, meta)

    align_germline(args, meta)
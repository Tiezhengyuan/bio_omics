'''
requirements:
    - table pdb_chain, meta_chain_faa
example: 
    - parallel -j4 python src/pipelines/imgt_msa.py ::: {0..46}
    - python src/pipelines/imgt_msa.py 46
functions:
    - collect data from tables and perform MSA
'''
from collections import defaultdict
from dotenv import load_dotenv
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import Dir, QueryComplex, AlignSeq, CalPhylo, \
    QueryComplex, ProcessSeq, ProcessPickle, BuildComplex

def antibody_seq(pdb_group, chain_type, cond):
    qc = QueryComplex(args['verbose'])
    query = f"""
        SELECT chain_id AS id, pdb_id, chain_seq As seq, pdb_file
        FROM view_antibody
        WHERE pdb_group={pdb_group}
            AND chain_seq IS NOT NULL
            AND chain_type IN {cond}
            AND CONCAT(chain_id, '{chain_type}') NOT IN (
                SELECT CONCAT(chain_id, chain_types)
                FROM msa_vregion
            );
    """
    rows = qc.list_data(query)
    return rows


def germline_human_vregion(pdb_group, args, meta):
    meta['pdb_group'] = int(pdb_group)
    print('Retrieve antibody + germline sequence...')
    pool = [
        ["H", "('H')"],
        ["K", "('K')"],
        ["L", "('L')"],
        ["KL", "('K', 'L')"],
    ]
    for chain_type, cond in pool:
        abs = antibody_seq(pdb_group, chain_type, cond)
        meta[chain_type] += len(abs)
        if abs:
            germlines = QueryComplex(args['verbose']).germline_seq(cond)
            for ab in abs:
                yield chain_type, ab, germlines

def align_germline(row_iter, args, meta):
    for chain_type, ab, germlines in row_iter:
        chain_id = ab['id']
        pdb_id = ab['pdb_id']
        outdir = os.path.join(args['outdir'], pdb_id[:2], pdb_id, 'MSA')
        Dir(outdir).init_dir()
        outprefix = os.path.join(outdir, f"{chain_id}.{chain_type}")

        # build fasta
        fa_file = outprefix + '.habs.faa'
        if args['overwrite'] or not os.path.isfile(fa_file):
            del ab['pdb_id']
            records = [ab,] + germlines
            ProcessSeq.to_fasta(records, fa_file)
            meta['msa_query_faa'] += 1
        else:
            meta['skip_msa_query_faa'] += 1
        
        # multiple sequence alignment
        aln_file = outprefix + '.habs.afa'
        if args['overwrite'] or not os.path.isfile(aln_file):
            AlignSeq(args['verbose']).muscle(fa_file, aln_file, large_msa=True)
            meta['muscle_alignment'] += 1
        else:
            meta['skip_muscle_alignment'] += 1

        pfile = outprefix + '.habs.pkl'
        if os.path.isfile(aln_file) and (args['overwrite'] or not os.path.isfile(pfile)):
            print('Try to calculate phylogentic distance...')
            data = CalPhylo().msa_distance(aln_file)
            if args['verbose']:
                print(data['distance'])
            if data:
                data.update({
                    'chain_type': chain_type,
                    'chain_id': chain_id,
                    'pdb_id': pdb_id,
                    'fa_file': fa_file,
                    'antibody_chain_pdb': ab['pdb_file'],
                })
                ProcessPickle(pfile).save_pickle(data)
                if args['verbose']:
                    print('Export phylogentic distiance to ', pfile)
                meta['calculate_msa_distance'] += 1
            else:
                meta['fail_msa_distance'] += 1
        else:
            meta['skip_msa_distance'] += 1
        
        record = (chain_id, chain_type, fa_file, aln_file, pfile)
        yield record
        

##############################
if __name__ == "__main__":
    meta = defaultdict(int)
    load_dotenv()
    args = {
        'verbose':True,
        'overwrite':False,
        'outdir': os.path.join(os.getenv("data_dir"), 'data'),
        'chunk_size': 10,
    }
    pdb_group = sys.argv[1]
    row_iter = germline_human_vregion(pdb_group, args, meta)

    print('Do MSA...')
    record_iter = align_germline(row_iter, args, meta)

    # put meta data to table msa_antibody
    bc = BuildComplex(record_iter, args['verbose'], args['chunk_size'])
    bc.insert_msa_vregion()

    print(meta)

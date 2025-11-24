'''
example: 
    - parallel -j4 python src/pipelines/export_chain_faa.py ::: {0..46}
    - python src/pipelines/export_chain_faa.py 46
functions:
    - retrieve pdb to fasta by chains
    - blastp against uniprot_sprot
'''
from dotenv import load_dotenv
from collections import defaultdict
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import Dir, Pdb, QueryComplex, BuildComplex, ProcessFasta

def retrieve_pdb(pdb_group, args, meta):
    # meta_pdb.pdb_file is absolute path
    query = f"""
        SELECT pdb_id, pdb_file
        FROM meta_pdb
        WHERE pdb_group = {pdb_group}
            AND pdb_id NOT IN (
                SELECT DISTINCT P.pdb_id
                FROM pdb_chain P
                LEFT JOIN meta_chain_faa M
                ON M.chain_id = P.chain_id
                WHERE M.faa_file IS NOT NULL
                    OR P.chain_seq is NULL
            );
    """
    rows = QueryComplex(args['verbose']).list_data(query)
    print(f"{len(rows)} PDB files ready...")
    meta['pdb_group'] = int(pdb_group)
    meta['rows'] += len(rows)
    return rows

def export_chain_faa(rows, args, meta):
    for row in rows:
        pdb_id = row['pdb_id']
        dir_endpoint = os.path.join('data', pdb_id[:2], pdb_id, 'chain_faa')
        outdir = os.path.join(args['data_dir'], dir_endpoint)
        Dir(outdir).init_dir()

        # load structure
        p = Pdb(row['pdb_file'], None, args['verbose'])
        p.load_structure()
        # retrieve chains
        chains = p.parse_chains()
        for chain_id, _pdb_id, model_no, chain_no, chain in chains:
            file_name = f"{chain_id}.faa"
            outfile = os.path.join(outdir, file_name)
            if args['overwrite'] or (not os.path.isfile(outfile)):
                chain_seq = p._chain_seq(chain)
                pf = ProcessFasta(outfile, args['verbose'])
                status = pf.one_seq(chain_seq, chain_id)
                meta[status] += 1
            else:
                meta['skip_faa'] += 1
            if os.path.isfile(outfile):
                relative_path = os.path.join(dir_endpoint, file_name)
                yield chain_id, relative_path
        if len(chains) == 0:
            meta['no_chain_seq'] += 1
    
if __name__ == "__main__":
    meta = defaultdict(int)
    load_dotenv()
    args = {
        'verbose': True,
        'overwrite': False,
        'data_dir': os.getenv('data_dir'),
        'chunk_size': 2,
    }
    pdb_group = sys.argv[1]
    rows = retrieve_pdb(pdb_group, args, meta)
    # export chain sequences in fasta to data/ in default
    record_iter = export_chain_faa(rows, args, meta)
    # insert
    bc = BuildComplex(record_iter, args['verbose'], args['chunk_size'])
    bc.insert_meta_chain_faa()

    print(meta)

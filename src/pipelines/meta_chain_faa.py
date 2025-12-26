'''
example: 
    - parallel -j4 python src/pipelines/meta_chain_faa.py ::: {0..116}
    - python src/pipelines/meta_chain_faa.py 0
functions:
    - retrieve pdb to fasta by chains
    - used by aligners, namely blastp
    - used by predictors, namely alphafold, etc
'''

from helper import *
from bioomics import Pdb, ProcessFasta

def retrieve_pdb(pdb_group, args, meta):
    table_name = args['table_name']
    query = f"""
        SELECT pdb_id, relative_path
        FROM meta_pdb
        WHERE pdb_group = {pdb_group}
            AND pdb_id NOT IN (
                SELECT DISTINCT pdb_id FROM {table_name}
            );
    """
    rows = QueryComplex(args['verbose']).list_data(query)
    print(f"{len(rows)} PDB files ready...")
    meta['pdb_group'] = int(pdb_group)
    meta['rows'] += len(rows)
    return rows

def export_chain_faa(rows, args, meta):
    for row in rows:
        # load structure
        pdb_id = row['pdb_id']
        pdb_file = os.path.join(args['pdb_dir'], row['relative_path'])
        p = Pdb(pdb_file, None, args['verbose'])
        p.load_structure()

        # retrieve chains
        chains = p.parse_chains()
        if chains:
            endpoint = os.path.join(args['table_name'], pdb_id[:2], pdb_id)
            outdir = os.path.join(args['data_dir'], endpoint)
            Dir(outdir).init_dir()
            for chain_id, _pdb_id, model_no, chain_no, chain in chains:
                file_name = f"{chain_id}.faa"
                outfile = os.path.join(outdir, file_name)
                if args['overwrite'] or (not os.path.isfile(outfile)):
                    chain_seq = p._chain_seq(chain)
                    if chain_seq:
                        pf = ProcessFasta(outfile, args['verbose'])
                        status = pf.one_seq(chain_seq, chain_id)
                        meta[status] += 1
                    else:
                        meta['empty'] += 1
                else:
                    meta['skip_faa'] += 1
                # build record
                if os.path.isfile(outfile):
                    relative_path = os.path.join(endpoint, file_name)
                    yield pdb_id, chain_id, relative_path
                else:
                    yield pdb_id, chain_id, None
        else:
            meta['no_chain_seq'] += 1
            yield pdb_id, None, None
    
if __name__ == "__main__":
    args = {
        'verbose': True,
        'overwrite': False,
        'pdb_dir': os.getenv('pdb_dir'),
        'data_dir': os.getenv('data_dir'),
        'chunk_size': 50,
        'table_name': 'meta_chain_faa',
    }
    pdb_group = sys.argv[1]
    rows = retrieve_pdb(pdb_group, args, meta)

    # export chain sequences in fasta to data/ in default
    record_iter = export_chain_faa(rows, args, meta)

    # insert
    bc = BuildComplex(record_iter, args['verbose'], args['chunk_size'])
    bc.insert_batch_records(args['table_name'])

    print(meta)

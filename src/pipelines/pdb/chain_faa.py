'''
example: 
    - parallel -j4 python app.py chain_faa ::: {0..240}
    - python app.py chain_faa 0
functions:
    - retrieve pdb to fasta by chains
    - used by aligners, namely blastp
    - used by predictors, namely alphafold, etc
'''

from src.helper import *
from bioomics import Pdb, ProcessFasta

@header
def retrieve_pdb(params, meta):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        SELECT pdb_id, relative_path
        FROM meta_pdb
        WHERE pdb_group = {pdb_group}
            AND pdb_id NOT IN (
                SELECT DISTINCT pdb_id FROM {table_name}
            );
    """
    rows = QueryComplex(params['verbose']).list_data(query)
    print(f"{len(rows)} PDB files ready...")
    meta['rows'] += len(rows)
    return rows

def export_chain_faa(rows, params, meta):
    for row in rows:
        # load structure
        pdb_id = row['pdb_id']
        pdb_file = os.path.join(params['rawdata_pdb_dir'], row['relative_path'])
        p = Pdb(pdb_file, None, params['verbose'])
        p.load_structure()

        # retrieve chains
        chains = p.parse_chains()
        if chains:
            endpoint = os.path.join(params['table_name'], pdb_id[:2], pdb_id)
            outdir = os.path.join(params['output_pdb_dir'], endpoint)
            Path(outdir).mkdir(parents=True, exist_ok=True)

            for item in chains:
                chain_id = item['chain_id']
                file_name = f"{chain_id}.faa"
                outfile = os.path.join(outdir, file_name)
                if params['overwrite'] or (not os.path.isfile(outfile)):
                    chain_seq = p._chain_seq(item['chain'])
                    if chain_seq:
                        pf = ProcessFasta(outfile, params['verbose'])
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
    params.update({
        'overwrite': False,
        'chunk_size': 50,
        'table_name': 'chain_faa',
    })
    
    rows = retrieve_pdb(params, meta)

    # export chain sequences in fasta to data/ in default
    record_iter = export_chain_faa(rows, params, meta)

    # insert
    bc = BuildComplex(params['verbose'], params['chunk_size'])
    bc.insert_batch_records(record_iter, params['table_name'])

    footer(meta)
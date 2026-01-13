'''
example: 
    - parallel -j4 python src/pipelines/meta_chain_pdb.py ::: {0..117}
    - python src/pipelines/meta_chain_pdb.py 0
functions:
    - retrieve pdb to pdb by chain
    - data are required by B-factor, DSSP analysis
note:
    - only include model.id=0 currently
'''

from helper import *
from bioomics import Pdb

def retrieve_pdb(pdb_group, args, meta):
    table_name = args['table_name']
    query = f"""
        SELECT pdb_id, relative_path
        FROM meta_pdb
        WHERE pdb_group = {pdb_group}
            AND pdb_id NOT IN (
                SELECT pdb_id FROM {table_name}
            );
    """
    rows = QueryComplex(args['verbose']).list_data(query)
    meta['pdb_group'] = int(pdb_group)
    meta['rows'] += len(rows)
    print(f"{len(rows)} PDB files ready for processing...")
    return rows

def pull_chain_pdb(rows, args, meta):
    for row in rows:
        pdb_id = row['pdb_id']
        pdb_file = os.path.join(args['pdb_dir'], row['relative_path'])
        try:
            p = Pdb(pdb_file, None, args['verbose'])
            p.load_structure()

            model_chains = p.parse_model_chains()
            model_no = 0
            if model_no in model_chains:
                endpoint = os.path.join(args['table_name'], pdb_id[:2], pdb_id)
                outdir = os.path.join(args['data_dir'], endpoint)
                Dir(outdir).init_dir()

                for item in model_chains[model_no]:
                    chain_id = item['chain_id']
                    file_name = f'{chain_id}.pdb'
                    outfile = os.path.join(outdir, file_name)
                    if args['overwrite'] or not os.path.isfile(outfile):
                        chain_no = item['chain_no']
                        p.to_chain_pdb(model_no, [chain_no,], outfile)
                        meta['chains'] += 1
                    else:
                        meta['skip'] += 1
                    
                    if os.path.isfile(outfile):
                        relative_path = os.path.join(endpoint, file_name)
                        yield pdb_id, chain_id, relative_path
                    else:
                        yield pdb_id, chain_id, None
            else:
                meta['null'] += 1
                yield pdb_id, '-', None
        except Exception as e:
            print(f"error={e}")
            meta['fail_parsing'] += 1
            yield pdb_id, None, None


if __name__ == "__main__":
    args = {
        'verbose': True,
        'overwrite': False,
        'chunk_size': 100,
        'pdb_dir': os.getenv("pdb_dir"),
        'data_dir': os.getenv("data_dir"),
        'table_name': 'meta_chain_pdb',
    }
    pdb_group = sys.argv[1]
    rows = retrieve_pdb(pdb_group, args, meta)
    record_iter = pull_chain_pdb(rows, args, meta)
    # insert
    bc = BuildComplex(record_iter, args['verbose'], args['chunk_size'])
    bc.insert_batch_records(args['table_name'])

    print(meta)

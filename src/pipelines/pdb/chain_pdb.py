'''
example: 
    - parallel -j8 python app.py chain_pdb ::: {0..240}
    - python app.py chain_pdb 0
functions:
    - retrieve pdb to pdb by chain
    - data are required by B-factor, DSSP analysis
note:
    - only include model.id=0 currently
'''

from src.helper import *
from bioomics import Pdb

@header
def retrieve_pdb(params, meta):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        SELECT pdb_id, relative_path
        FROM meta_pdb
        WHERE pdb_group = {pdb_group}
            AND pdb_id NOT IN (
                SELECT pdb_id FROM {table_name}
            );
    """
    rows = QueryComplex(params['verbose']).list_data(query)
    meta['pdb_group'] = int(pdb_group)
    meta['rows'] += len(rows)
    print(f"{len(rows)} PDB files ready for processing...")
    return rows

def pull_chain_pdb(rows, params, meta):
    for row in rows:
        pdb_id = row['pdb_id']
        pdb_file = os.path.join(params['rawdata_pdb_dir'], row['relative_path'])
        p = Pdb(pdb_file, None, params['verbose'])
        p.load_structure()

        model_chains = p.parse_model_chains()
        model_no = 0
        if model_no in model_chains:
            endpoint = os.path.join(params['table_name'], pdb_id[:2], pdb_id)
            outdir = os.path.join(params['output_pdb_dir'], endpoint)
            Path(outdir).mkdir(parents=True, exist_ok=True)

            for item in model_chains[model_no]:
                chain_id = item['chain_id']
                file_name = f'{chain_id}.pdb'
                outfile = os.path.join(outdir, file_name)
                if params['overwrite'] or not os.path.isfile(outfile):
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

if __name__ == "__main__":
    params.update({
        'overwrite': False,
        'chunk_size': 100,
        'table_name': 'chain_pdb',
    })
    
    rows = retrieve_pdb(params, meta)
    record_iter = pull_chain_pdb(rows, params, meta)
    # insert
    bc = BuildComplex(params['verbose'], params['chunk_size'])
    bc.insert_batch_records(record_iter, params['table_name'])

    print("Try to validate files...")
    validate_data('chain_pdb', params, meta)
    footer(meta)

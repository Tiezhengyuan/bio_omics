'''
example: 
    - parallel -j8 bomics chain_pdb ::: {0..241}
    - bomics chain_pdb 0
functions:
    - retrieve pdb to pdb by chain
    - data are required by B-factor, DSSP analysis
note:
    - only include model.id=0 currently
'''

from src.helper import *
from bioomics import Pdb

@header
def retrieve(params, meta):
    table_name = params['table_name']
    pdb_group = params['pdb_group']
    query = f"""
        SELECT c.pdb_id, c.chain_id, m.relative_pdb
        FROM pdb_chainid   c
        LEFT JOIN meta_pdb m ON c.pdb_id = m.pdb_id
        WHERE m.pdb_group = {pdb_group}
            AND c.chain_id NOT IN (
                SELECT chain_id FROM {table_name}
            )
    """
    rows = QueryComplex(params['verbose']).list_data(query)
    meta['rows'] = len(rows)
    if len(rows) > 0:
        df = pd.DataFrame(rows)
        g = df.groupby(['pdb_id', 'relative_pdb'])
        for (pdb_id, relative_pdb), sub in g:
            pdb_file = params['rawdata_dir'] / relative_pdb
            if pdb_file.is_file():
                endpoint = Path(params['table_name']) / pdb_id[:2] / pdb_id
                outdir = params['output_pdb_dir'] / endpoint
                outdir.mkdir(parents=True, exist_ok=True)
                chain_ids = sub['chain_id'].to_list()
                yield pdb_id, pdb_file, chain_ids, endpoint, outdir

def build(row_iter, params, meta):
    for pdb_id, pdb_file, chain_ids, endpoint, outdir in row_iter:
        p = Pdb(pdb_file, None, params['verbose'])
        p.load_structure()
                
        tag = 0 
        for item in p.parse_chains():
            chain_id = item['chain_id']
            if chain_id in chain_ids:
                meta['records'] += 1
                file_name = f'{chain_id}.pdb'
                outfile = outdir / file_name
                if params['overwrite'] or not outfile.is_file():
                    chain_no = item['chain_no']
                    p.to_chain_pdb(item['model_no'], [chain_no,], outfile)
                # build record
                relative_pdb = str(endpoint / file_name) if outfile.is_file() else None
                yield (pdb_id, chain_id, relative_pdb)
            else:
                meta['omit_chains'] += 1
            tag = 1
        if tag == 0:
            if params['verbose']:
                print(f"Warning: no chains {pdb_id} are detectin in {pdb_file}")
            yield (pdb_id, None, None)
            meta['null'] += 1
            

if __name__ == "__main__":
    params.update({
        'overwrite': False,
        'chunk_size': 100,
        'table_name': 'chain_pdb',
    })
    print("Try to validate files...")
    validate_data('relative_pdb', params, meta)
    
    rows = retrieve(params, meta)
    record_iter = build(rows, params, meta)
    # insert
    bc = BuildComplex(params['verbose'], params['chunk_size'])
    bc.insert_batch_records(record_iter, params['table_name'])

    footer(meta)

'''
requirements:
    - The table pdb_chainid is ready.
example:
    - parallel -j8 bomics pdb_chain_seq ::: {0..241}
    - bomics pdb_chain_seq 0
function:
    - assign unique id to chain protein sequences
    - export to <chain_id>.faa
'''
import pandas as pd
from src.helper import *
from bioomics import Pdb, ProcessFasta

@header
def retrieve(params, meta):
    table_name = params['table_name']
    pdb_group = params['pdb_group']
    query = f"""
        SELECT c.pdb_id, c.chain_id, m.relative_pdb
        FROM pdb_chainid   c
        LEFT JOIN meta_pdb m ON c.pdb_id = m.pdb_id
        WHERE m.pdb_group = {pdb_group}
            AND NOT EXISTS (
                SELECT 1 FROM {table_name} a
                WHERE c.chain_id = a.chain_id
            )
    """
    rows = QueryComplex(params['verbose']).list_data(query)
    meta['rows'] = len(rows)
    print(f"Detect chain ids: {len(rows)}")
    if len(rows) > 0:
        df = pd.DataFrame(rows)
        g = df.groupby(['pdb_id', 'relative_pdb'])
        for (pdb_id, relative_pdb), sub in g:
            pdb_file = params['rawdata_dir'] / relative_pdb
            if pdb_file.is_file():
                endpoint = Path(params['table_name']) / pdb_id[:2] / pdb_id
                outdir = params['output_pdb_dir'] / endpoint
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
                # the column pro_seq is not null
                chain_seq = p._chain_seq(item['chain'])
                if chain_id and len(chain_seq) > 0:
                    tag = 1
                    # save seq as fasta
                    outdir.mkdir(parents=True, exist_ok=True)
                    file_name = f"{chain_id}.faa"
                    fa_file = outdir / file_name
                    if params['overwrite'] or not fa_file.is_file():
                        ProcessFasta(fa_file).one_seq(chain_seq, chain_id)
                    # build record
                    relative_faa = str(endpoint / file_name) if fa_file.is_file() else None
                    yield (pdb_id, chain_id, chain_seq, relative_faa)
                else:
                    yield (pdb_id, chain_id, None, None)
            else:
                meta['omit_chains'] += 1
        if tag == 0:
            if params['verbose']:
                print(f"Warning: no chains {pdb_id} are detectin in {pdb_file}")
            yield (pdb_id, None, None, None)
            meta['null'] += 1

##################################
if __name__ == "__main__":
    params.update({
        'overwrite': False,
        'table_name': 'pdb_chain_seq',
    })
    validate_data('relative_faa', params, meta)
    
    row_iter = retrieve(params, meta)
    record_iter = build(row_iter, params, meta)

    print(f"Try to insert PDB data into table {params['table_name']}...")
    bc = BuildComplex(params['verbose'])
    bc.insert_batch_records(record_iter, params['table_name'])
   
    footer(meta)
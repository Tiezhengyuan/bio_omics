'''
requirments:
    - chain_pdb and pdb_chain_seq
    - label_chain
example: 
    - parallel -j8 python app.py equal_seq_rmsd ::: {0..240}
    - python app.py equal_seq_rmsd 0
functions:
    - calculate distance of two pdbs of which sequences are identical
    - put data to table equal_seq_rmsd
'''
from src.helper import *
from bioomics import EvalPred

def pull_chain_pdb(params, meta):
    query = f"""
        SELECT chain_id, chain_pdb FROM chain_pdb
        WHERE chain_pdb IS NOT NULL
    """
    rows = QueryComplex(params['verbose']).list_data(query)
    meta['num_chains'] = len(rows)
    df = pd.DataFrame(rows)
    res = dict(zip(df['chain_id'], df['chain_pdb']))
    return res

def pull_equal_rmsd(params, meta):
    table_name = params['table_name']
    query = f"""
        SELECT chain_id, equal_chain_id
        FROM {table_name}
    """
    rows = QueryComplex(params['verbose']).list_data(query)
    pairs = [set((i['chain_id'], i['equal_chain_id'])) for i in rows]
    meta['equal_pairs'] += len(rows)
    return pairs

def pull_equal_ids(pairs, params, meta):
    pdb_group = params['pdb_group']
    query = f"""
        WITH tmp AS (
            SELECT first_chain_id 
            FROM label_chain   c
            LEFT JOIN meta_pdb p ON c.pdb_id = p.pdb_id
            WHERE p.pdb_group = {pdb_group}    
        )
        SELECT chain_id, first_chain_id
        FROM pdb_chain_seq
        WHERE chain_seq IS NOT NULL
            AND first_chain_id IN (
                SELECT first_chain_id FROM tmp
            )
    """
    rows = QueryComplex(params['verbose']).list_data(query)
    if rows:
        df = pd.DataFrame(rows)
        for _, sub in df.groupby('first_chain_id'):
            if len(sub) > 1:
                meta['num_equal_seq'] += 1
                ids = sub['chain_id'].tolist()
                combo_iter = itertools.combinations(ids, 2)
                # exclude records already in table equal_seq_rmsd
                for (a, b) in combo_iter:
                    if set((a, b)) not in pairs:
                        yield (a, b)

def calculate_rmsd(ids_pool, chain_pdb, args, meta):
    ep = EvalPred()
    for (id1, id2) in ids_pool:
        relative1 = chain_pdb.get(id1)
        relative2 = chain_pdb.get(id2)
        rmsd = None
        if relative1 and relative2:
            pdb1 = os.path.join(args['output_pdb_dir'], relative1)
            pdb2 = os.path.join(args['output_pdb_dir'], relative2)
            if os.path.isfile(pdb1) and os.path.isfile(pdb2):
                rmsd = ep.rmsd(pdb1, pdb2)
                if rmsd:
                    meta['succeed'] += 1
                else:
                    meta['fail'] += 1
        yield (id1, id2, rmsd, relative1, relative2)

if __name__ == "__main__":
    params.update({
        'chunk_size': 20,
        'table_name': 'equal_seq_rmsd',
    })

    # prepare data
    chain_pdb = pull_chain_pdb(params, meta)
    pairs = pull_equal_rmsd(params, meta)
    ids_pool = pull_equal_ids(pairs, params, meta)
    
    # calculate atom distance
    record_iter = calculate_rmsd(ids_pool, chain_pdb, params, meta)

    # put meta data to table
    bc = BuildComplex(params['verbose'], params['chunk_size'])
    bc.insert_batch_records(record_iter, params['table_name'])

    footer(meta)
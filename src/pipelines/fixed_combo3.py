'''
requirements:
    - The table fixed_pdb is ready.
example:
    - parallel -j2 python src/pipelines/fixed_combo3.py ::: {0..116}
    - python src/pipelines/fixed_combo3.py 0
function:
    - split fixed PDB of complex to chain PDB.
note: 
    - only consider model.id = 0
'''

from itertools import combinations
from helper import *

def retrieve_data(pdb_group, args, meta):
    table_name = args['table_name']
    query = f"""
        SELECT F.pdb_id, F.relative_path
        FROM fixed_pdb F
        LEFT JOIN meta_pdb P ON F.pdb_id = P.pdb_id
        WHERE P.pdb_group = {pdb_group}
            AND F.pdb_id NOT IN (
                SELECT pdb_id FROM {table_name}
            )
            AND F.relative_path IS NOT NULL
    ;"""
    rows = QueryComplex(args['verbose']).list_data(query)
    meta['pdb_group'] = int(pdb_group)
    meta['rows'] = len(rows)
    return rows

def pull_chain_pdb(rows, args, meta):
    combo = int(args['combo'])
    for row in rows:
        pdb_file = os.path.join(args['pdb_dir'], row['relative_path'])
        p = Pdb(pdb_file, None, args['verbose'])
        p.load_structure()

        pdb_id = row['pdb_id']
        model_no = 0
        rec0 = (pdb_id, model_no,)
        chains = p.model_chains_seq(model_no)
        chain_no_pool = list(chains)
        if len(chain_no_pool) >= combo:
            endpoint = os.path.join(args['table_name'], pdb_id[:2], pdb_id) 
            outdir = os.path.join(args['pdb_dir'], endpoint)
            Dir(outdir).init_dir()

            chain_no_pool.sort()
            for _ids in combinations(chain_no_pool, combo):
                chain_combo = '-'.join(list(_ids))
                combo_id = f"{pdb_id}_{model_no}_{chain_combo}"
                file_name = f"{combo_id}.pdb"
                outfile = os.path.join(outdir, file_name)
                if args['overwrite'] or not os.path.isfile(outfile):
                    status = p.to_chain_pdb(model_no, list(_ids), outfile)
                    meta[status] += 1
                else:
                    meta['skip'] += 1
                # prepare record
                relative_path = os.path.join(endpoint, file_name) if os.path.isfile(outfile) else None
                yield rec0 + (combo_id, chain_combo, relative_path)
        else:
            yield rec0 + (None,) * 3

##################################
if __name__ == "__main__":
    pdb_group = sys.argv[1]
    args = {
        'verbose': True,
        'overwrite': False,
        'chunk_size': 50,
        'pdb_dir': os.getenv('pdb_dir'),
        'table_name': 'fixed_combo3',
        'combo': 3,
    }
    rows = retrieve_data(pdb_group, args, meta)

    if len(rows) > 0:
        print(f"Try to insert {len(rows)} data into {args['table_name']}...")
        record_iter = pull_chain_pdb(rows, args, meta)
        bc = BuildComplex(record_iter, args['verbose'], args['chunk_size'])
        bc.insert_batch_records(args['table_name'])

    print(meta)


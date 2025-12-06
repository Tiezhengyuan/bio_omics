'''
requirements:
    - The table fixed_pdb is ready.
example:
    - parallel -j2 python src/pipelines/fixed_combo3.py ::: {0..116}
    - python src/pipelines/fixed_combo3.py 0
function:
    - split fixed PDB of complex to chain PDB.
'''
from itertools import combinations
from dotenv import load_dotenv
from collections import defaultdict
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import QueryComplex, BuildComplex, Dir, Pdb


def retrieve_data(pdb_group, args, meta):
    table_name = args['table_name']
    query = f"""
        SELECT F.pdb_id, F.pdb_file
        FROM fixed_pdb F
        LEFT JOIN meta_pdb P
        ON F.pdb_id = P.pdb_id
        WHERE P.pdb_group = {pdb_group}
            AND F.pdb_id NOT IN (
                SELECT pdb_id FROM {table_name}
            )
            AND F.pdb_file IS NOT NULL
    ;"""
    rows = QueryComplex(args['verbose']).list_data(query)
    meta['pdb_group'] = int(pdb_group)
    meta['rows'] = len(rows)
    return rows

def pull_chain_pdb(rows, args, meta):
    combo = int(args['combo'])
    for row in rows:
        pdb_id = row['pdb_id']
        pdb_file = os.path.join(args['pdb_dir'], row['pdb_file'])
        p = Pdb(pdb_file, None, args['verbose'])

        chain_no_pool = None
        try:
            p.load_structure()
            chain_no_pool = p.get_chain_no()
        except Exception  as e:
            print(e)
            continue

        if chain_no_pool:
            dir_endpoint = os.path.join(args['table_name'], pdb_id[:2], pdb_id) 
            outdir = os.path.join(args['pdb_dir'], dir_endpoint)
            Dir(outdir).init_dir()
            for model_no, _pool in chain_no_pool.items():
                _pool.sort()
                if len(_pool) >= combo:
                    for _ids in combinations(_pool, combo):
                        chain_combo = '-'.join(list(_ids))
                        file_name = f"{pdb_id}_{model_no}_{chain_combo}.pdb"
                        outfile = os.path.join(outdir, file_name)
                        if args['overwrite'] or not os.path.isfile(outfile):
                            status = p.to_chain_pdb(list(_ids), outfile)
                            meta[status] += 1
                        else:
                            meta['skip'] += 1
                        # prepare record
                        if os.path.isfile(outfile):
                            relative_path = os.path.join(dir_endpoint, file_name)
                            yield pdb_id, model_no, chain_combo, relative_path
                        else:
                            yield pdb_id, model_no, chain_combo, None
                else:
                    yield pdb_id, model_no, None, None
        else:
            yield pdb_id, None, None, None

##################################
if __name__ == "__main__":
    pdb_group = sys.argv[1]
    load_dotenv()
    args = {
        'verbose': True,
        'overwrite': True,
        'pdb_dir': os.getenv('pdb_dir'),
        'chunk_size': 50,
        'combo': 3,
        'table_name': 'fixed_combo3',
    }
    meta = defaultdict(int)
    rows = retrieve_data(pdb_group, args, meta)

    if len(rows) > 0:
        print(f"Try to insert {len(rows)} data into {args['table_name']}...")
        record_iter = pull_chain_pdb(rows, args, meta)
        bc = BuildComplex(record_iter, args['verbose'], args['chunk_size'])
        bc.insert_batch_records(args['table_name'])

    print(meta)


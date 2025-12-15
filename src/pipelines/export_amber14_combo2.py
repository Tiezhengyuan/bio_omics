'''
requirements:
    - The table relaxed_pdb_amber14 is ready.
example:
    - parallel -j2 python src/pipelines/export_amber14_combo2.py ::: {0..116}
    - python src/pipelines/export_amber14_combo2.py 0
function:
    - split relaxed PDB of complex to chain PDB.
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
        FROM relaxed_pdb_amber14 F
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

def build_combo(rows, args, meta):
    combo = int(args['combo'])
    for row in rows:
        pdb_id = row['pdb_id']
        pdb_file = os.path.join(args['pdb_dir'], row['pdb_file'])
        p = Pdb(pdb_file, None, args['verbose'])

        chains_seq = None
        try:
            p.load_structure()
            chains_seq = p.model_chains_seq()
        except Exception  as e:
            print(e)
            continue

        if chains_seq:
            dir_endpoint = os.path.join(args['table_name'], pdb_id[:2], pdb_id) 
            outdir = os.path.join(args['data_dir'], dir_endpoint)
            Dir(outdir).init_dir()
            for model_no, seq_pool in chains_seq.items():
                if len(seq_pool) >= combo:
                    _pool = list(seq_pool)
                    _pool.sort()
                    for _ids in combinations(_pool, combo):
                        chain_combo = '-'.join(list(_ids))
                        combo_id = f"{pdb_id}_{model_no}_{chain_combo}"
                        # export data
                        out_pdb = os.path.join(outdir, f"{combo_id}.pdb")
                        out_csv = os.path.join(outdir, f"{combo_id}.csv")
                        if args['overwrite'] or not os.path.isfile(out_pdb):
                            # to *.faa
                            combo_seq = ':'.join([seq_pool[i] for i in _ids])
                            text = ['id,sequence', combo_id + ',' + combo_seq]
                            with open(out_csv, 'w') as f:
                                f.write('\n'.join(text))
                            # to *.pdb
                            status = p.to_chain_pdb(list(_ids), out_pdb)
                            meta[status] += 1
                        else:
                            meta['skip'] += 1
                        # prepare record
                        if os.path.isfile(out_pdb):
                            endpoint_prefix = os.path.join(dir_endpoint, combo_id)
                            yield pdb_id, combo_id, model_no, chain_combo, endpoint_prefix
                        else:
                            yield pdb_id, combo_id, model_no, chain_combo, None
                else:
                    yield pdb_id, None, model_no, None, None
        else:
            yield pdb_id, None, None, None, None

##################################
if __name__ == "__main__":
    pdb_group = sys.argv[1]
    load_dotenv()
    args = {
        'verbose': True,
        'overwrite': True,
        'pdb_dir': os.getenv('pdb_dir'),
        'data_dir': os.getenv('data_dir'),
        'chunk_size': 50,
        'combo': 2,
        'table_name': 'relaxed_amber14_combo2',
    }
    meta = defaultdict(int)
    rows = retrieve_data(pdb_group, args, meta)

    if len(rows) > 0:
        print(f"Try to insert {len(rows)} data into {args['table_name']}...")
        record_iter = build_combo(rows, args, meta)
        bc = BuildComplex(record_iter, args['verbose'], args['chunk_size'])
        bc.insert_batch_records(args['table_name'])

    print(meta)


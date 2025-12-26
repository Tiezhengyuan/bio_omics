'''
requirements:
    - The table meta_pdb is ready.
example:
    - parallel -j4 python src/pipelines/meta_combo2_seq.py ::: {0..116}
    - python src/pipelines/meta_combo2_seq.py 0
function:
    - split PDB of complex to paired chains PDB.
note:
    - only consider model.id == 0
'''
from itertools import combinations
from helper import *
from bioomics import ProcessSeq

def retrieve_data(pdb_group, args, meta):
    table_name = args['table_name']
    query = f"""
        SELECT pdb_id, relative_path
        FROM meta_pdb
        WHERE pdb_group = {pdb_group}
            AND pdb_id NOT IN (
                SELECT pdb_id FROM {table_name}
            )
            AND relative_path IS NOT NULL
    ;"""
    rows = QueryComplex(args['verbose']).list_data(query)
    meta['pdb_group'] = int(pdb_group)
    meta['rows'] = len(rows)
    return rows

def pull_chain_pdb(rows, args, meta):
    combo = int(args['combo'])
    for row in rows:
        pdb_id = row['pdb_id']
        pdb_file = os.path.join(args['pdb_dir'], row['relative_path'])
        try:
            p = Pdb(pdb_file, None, args['verbose'])
            p.load_structure()

            # only consider model 0
            model_no = 0
            chains_seq = p.model_chains_seq(model_no)
            ids = list(chains_seq)
            ids.sort()
            if len(ids) >= combo:
                endpoint = os.path.join(args['table_name'], pdb_id[:2], pdb_id) 
                outdir = os.path.join(args['data_dir'], endpoint)
                Dir(outdir).init_dir()

                for _ids in combinations(ids, combo):
                    chain_combo = '-'.join(list(_ids))
                    combo_id = f"{pdb_id}_{model_no}_{chain_combo}"
                    file_name = f"{combo_id}.csv"
                    outfile = os.path.join(outdir, file_name)
                    if args['overwrite'] or not os.path.isfile(outfile):
                        combo_seq = {
                            combo_id: [chains_seq[i] for i in _ids],
                        }
                        status = ProcessSeq.to_colabfold_csv(combo_seq, outfile)
                        meta[status] += 1
                    else:
                        meta['skip'] += 1
                    # prepare record
                    rec = [pdb_id, model_no, combo_id, chain_combo, None]
                    if os.path.isfile(outfile):
                        rec[-1] = os.path.join(endpoint, file_name)
                    yield tuple(rec)
            else:
                # PDB parsing is good, but no model 0 or chains < combo
                yield (pdb_id, model_no) + (None,) * 3
        except Exception  as e:
            print(f"Error: Can't retrieve PDB from {pdb_file}, error={e}")
            yield (pdb_id,) + (None,) * 4

##################################
if __name__ == "__main__":
    pdb_group = sys.argv[1]
    args = {
        'verbose': True,
        'overwrite': True,
        'pdb_dir': os.getenv('pdb_dir'),
        'data_dir': os.getenv('data_dir'),
        'chunk_size': 50,
        'combo': 2,
        'table_name': 'meta_combo2_seq',
    }
    rows = retrieve_data(pdb_group, args, meta)

    if len(rows) > 0:
        print(f"Try to insert {len(rows)} data into {args['table_name']}...")
        record_iter = pull_chain_pdb(rows, args, meta)
        bc = BuildComplex(record_iter, args['verbose'], args['chunk_size'])
        bc.insert_batch_records(args['table_name'])

    print(meta)


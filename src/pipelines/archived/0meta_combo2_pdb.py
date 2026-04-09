'''
requirements:
    - The table meta_pdb is ready.
    - meta_chain_pdb is ready
example:
    - parallel -j4 python src/pipelines/meta_combo2_pdb.py ::: {0..117}
    - python src/pipelines/meta_combo2_pdb.py 0
function:
    - split PDB of complex to paired chains PDB.
note:
    - only consider model.id == 0
'''
from itertools import combinations
from helper import *

def retrieve_data(pdb_group, args, meta):
    table_name = args['table_name']
    query = f"""
        SELECT C.pdb_id, I.model_no, I.chain_no, C.chain_id, C.chain_pdb, P.relative_path
        FROM meta_chain_pdb C
        LEFT JOIN pdb_chainid   I ON C.chain_id=I.chain_id
        LEFT JOIN meta_pdb      P ON C.pdb_id=P.pdb_id
        WHERE P.pdb_group = {pdb_group}
            AND P.relative_path IS NOT NULL
            AND C.model_no = 0
            AND C.chain_id IS NOT NULL
            AND C.pdb_id NOT IN (
                SELECT pdb_id FROM {table_name}
            )
    ;"""
    df = QueryComplex(args['verbose']).list_data(query, True)
    meta['pdb_group'] = int(pdb_group)
    meta['rows'] = len(df)
    print(f"{len(df)} PDB files ready for data processing...")

    if len(df) > 0:
        g = df.groupby(['pdb_id', 'relative_path', 'model_no'])
        for (pdb_id, relative_path, model_no), sub in g:
            if len(sub) >= args['combo']:
                pdb_file = os.path.join(args['pdb_dir'], relative_path)
                if os.path.isfile(pdb_file):
                    yield {
                        'input_pdb': pdb_file,
                        'pdb_id': pdb_id,
                        'model_no': model_no,
                        'chain_pdbs': sub,
                    }

def pull_chain_pdb(input, args, meta):
    combo = int(args['combo'])
    pdb_id = input['pdb_id']
    endpoint = os.path.join(args['table_name'], pdb_id[:2], pdb_id) 
    outdir = os.path.join(args['data_dir'], endpoint)
    Dir(outdir).init_dir()

    try:
        p = Pdb(input['input_pdb'], None, args['verbose'])
        p.load_structure()

        model_no = input['model_no']
        ids = input['chain_pdbs']['chain_no']
        for _ids in combinations(ids, combo):
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
            rec = [pdb_id, model_no, combo_id, chain_combo, None, None]
            if os.path.isfile(outfile):
                rec[-1] = os.path.join(endpoint, file_name)
            chain_pdbs = []
            for chain_no in ids:
                chain_id = f"{pdb_id}_{model_no}_{chain_no}_0"
                chain_pdb = os.path.

            yield tuple(rec)

    except Exception  as e:
        print(f"Error: Can't retrieve PDB from {pdb_file}, error={e}")
        yield (pdb_id,) + (None,) * 5

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
        'table_name': 'meta_combo2_pdb',
    }
    input_iter = retrieve_data(pdb_group, args, meta)

    print(f"Try to insert {len(rows)} data into {args['table_name']}...")
    record_iter = pull_chain_pdb(input_iter, args, meta)
    bc = BuildComplex(record_iter, args['verbose'], args['chunk_size'])
    bc.insert_batch_records(args['table_name'])

    print(meta)


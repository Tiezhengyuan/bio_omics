'''
requirements:
    - The table fixed_pdb is ready.
example:
    - parallel -j4 python src/pipelines/fixed_chain.py ::: {0..116}
    - python src/pipelines/fixed_chain.py 0
function:
    - split fixed PDB of complex to chain PDB.
'''

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
    for row in rows:
        pdb_id = row['pdb_id']
        dir_endpoint = os.path.join(args['table_name'], pdb_id[:2], pdb_id) 
        outdir = os.path.join(args['pdb_dir'], dir_endpoint)
        Dir(outdir).init_dir()

        pdb_file = os.path.join(args['pdb_dir'], row['relative_path'])
        p = Pdb(pdb_file, None, args['verbose'])
        p.load_structure()
        chains = p.parse_chains(pdb_id)
        if chains:
            for item in chains:
                chain_no = item['chain_no']
                model_no = item['model_no']
                chain_id = item['chain_id']
                file_name = f"{chain_id}.pdb"
                outfile = os.path.join(outdir, file_name)
                if args['overwrite'] or (not os.path.isfile(outfile)):
                    status = p.to_chain_pdb(model_no, [chain_no,], outfile)
                    meta[status] += 1
                else:
                    meta['skip'] += 1
                rec0 = (pdb_id, chain_id, model_no, chain_no)
                if os.path.isfile(outfile):
                    relative_path = os.path.join(dir_endpoint, file_name)
                    yield rec0 + (relative_path,)
                else:
                    yield rec0 + (None,)
        else:
            yield pdb_id, None, None, None, None

##################################
if __name__ == "__main__":
    args = {
        'verbose': True,
        'overwrite': False,
        'chunk_size': 50,
        'pdb_dir': os.getenv('pdb_dir'),
        'table_name': 'fixed_chain',
    }
    pdb_group = sys.argv[1]
    rows = retrieve_data(pdb_group, args, meta)

    if len(rows) > 0:
        print(f"Try to insert {len(rows)} data into {args['table_name']}...")
        record_iter = pull_chain_pdb(rows, args, meta)
        bc = BuildComplex(record_iter, args['verbose'], args['chunk_size'])
        bc.insert_batch_records(args['table_name'])

    print(meta)
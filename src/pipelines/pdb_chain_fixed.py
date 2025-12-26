#!/usr/bin/python
'''
requirements:
    - table meta_chain_pdb
example:
    - parallel -j4 python src/pipelines/pdb_chain_fixed.py ::: {0..116}
    - python src/pipelines/pdb_chain_fixed.py 0
function:
    - fix pdb files
    - insert data into table pdb_chain_fixed
note:
    - only consider model.id=0 and antibody for saving time.
'''

from helper import *
from bioomics import OpenmmPdbfixer

def retrieve_pdb(args, meta):
    pdb_group = args['pdb_group']
    table_name = args['table_name']
    query = f"""
        SELECT M.pdb_id, M.chain_id, M.relative_path
        FROM meta_chain_pdb M
        LEFT JOIN pdb_chainid C ON M.chain_id = C.chain_id
        LEFT JOIN meta_pdb P ON M.pdb_id = P.pdb_id
        WHERE P.pdb_group = {pdb_group}
            AND C.model_no = 0
            AND M.relative_path IS NOT NULL
            AND M.chain_id NOT IN (
                SELECT chain_id FROM {table_name}
            )
            AND M.pdb_id IN (
                SELECT pdb_id FROM view_antibody
            )
    ;"""
    rows = QueryComplex(args['verbose']).list_data(query)
    print(f"{len(rows)} PDB files ready...")
    meta['pdb_group'] = int(pdb_group)
    meta['rows'] += len(rows)
    return rows

def run_fixer(rows, args, meta):
    for row in rows:
        time0 = datetime.now()
        pdb_id = row['pdb_id']
        endpoint = os.path.join(args['table_name'], pdb_id[:2], pdb_id)
        outdir = os.path.join(args['pdb_dir'], endpoint)
        Dir(outdir).init_dir()
        
        chain_id = row['chain_id']
        file_name = chain_id + '.pdb'
        fixed_pdb_file = os.path.join(outdir, file_name)
        if args['overwrite'] or not os.path.isfile(fixed_pdb_file):
            # run pdbfixer
            input_pdb = os.path.join(args['data_dir'], row['relative_path'])
            op = OpenmmPdbfixer(input_pdb, args['verbose'])
            output_prefix = os.path.join(outdir, chain_id)
            info = op.fast_fixer(output_prefix)
            meta[info['status']] += 1
        else:
            meta['skip'] += 1
        
        # build record
        if os.path.isfile(fixed_pdb_file):
            meta['insert_fixed'] += 1
            relative_path = os.path.join(endpoint, file_name)
            yield chain_id, relative_path, time0
        else:
            meta['insert_null'] += 1
            yield chain_id, None, time0


if __name__ == "__main__":
    args = {
        'verbose': True,
        'overwrite': False,
        'chunk_size': 5,
        'pdb_dir': os.getenv('pdb_dir'),
        'data_dir': os.getenv('data_dir'),
        'table_name': 'pdb_chain_fixed',
    }
    args['pdb_group'] = sys.argv[1]
    rows = retrieve_pdb(args, meta)

    if len(rows) > 0:
        # fix pdb
        record_iter = run_fixer(rows, args, meta)
        # insertion
        bc = BuildComplex(record_iter, args['verbose'], args['chunk_size'])
        bc.insert_batch_records(args['table_name'])

    print(meta)
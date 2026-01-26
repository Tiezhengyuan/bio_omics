#!/usr/bin/python
'''
requirements:
    - table meta_chain_pdb
example:
    - parallel -j8 python src/pipelines/pdb_chain_fix.py ::: {0..116}
    - python src/pipelines/pdb_chain_fix.py 0
function:
    - fix pdb files by chain
    - insert data into table chain_fix
'''

from helper import *
from bioomics import OpenmmPdbfixer

def retrieve_pdb(args, meta):
    pdb_group = args['pdb_group']
    table_name = args['table_name']
    # meta_chain_pdb.pdb_file is absolute path
    query = f"""
        SELECT P.pdb_id, M.chain_id, M.pdb_file
        FROM meta_chain_pdb M
        LEFT JOIN pdb_chainid C
        ON M.chain_id = C.chain_id
        LEFT JOIN meta_pdb P
        ON C.pdb_id = P.pdb_id
        WHERE P.pdb_group = {pdb_group}
            AND M.chain_id NOT IN (
                SELECT chain_id FROM {table_name}
            );
    """
    qc = QueryComplex(args['verbose'])
    rows = qc.list_data(query)
    print(f"{len(rows)} PDB files ready...")
    meta['pdb_group'] = int(pdb_group)
    meta['rows'] += len(rows)
    return rows

def run_fixer(rows, args, meta):
    for row in rows:
        pdb_id = row['pdb_id']
        path_endpoint = os.path.join('data', pdb_id[:2], pdb_id, 'chain_fixed')
        outdir = os.path.join(args['data_dir'], path_endpoint)
        Dir(outdir).init_dir()

        chain_id = row['chain_id']
        output_prefix = os.path.join(outdir, chain_id)
        fixed_pdb_file = output_prefix + '.pdb'
        if args['overwrite'] or not os.path.isfile(fixed_pdb_file):
            # run pdbfixer
            input_pdb = os.path.join(args['data_dir'], row['pdb_file'])
            op = OpenmmPdbfixer(input_pdb, args['verbose'])
            info = op.fast_fixer(output_prefix)
            meta[info["status"]] += 1
        else:
            meta['skip'] += 1
        
        # build record
        if os.path.isfile(fixed_pdb_file):
            relative_path = os.path.join(path_endpoint, f"{chain_id}.pdb")
            yield chain_id, relative_path
        else:
            meta['null'] += 1
            yield chain_id, None

if __name__ == "__main__":
    args = {
        'verbose': True,
        'overwrite': True,
        'data_dir': os.getenv('data_dir'),
        'chunk_size': 5,
        'table_name': 'chain_fix',
    }
    args['pdb_group'] = sys.argv[1]
    rows = retrieve_pdb(args, meta)

    if len(rows) > 0:
        # fix pdb
        record_iter = run_fixer(rows, args, meta)
        # insertion
        bc = BuildComplex(record_iter, args['verbose'], args['chunk_size'])
        bc.insert_chain_fix(args['table_name'])

    print(meta)
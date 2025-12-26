#!/usr/bin/python
'''
requirements:
    - PDB raw data
    - table meta_pdb
example:
    - parallel -j4 python src/pipelines/fixed_pdb.py ::: {0..116}
    - python src/pipelines/fixed_pdb.py 0
function:
    - fix pdb files
    - insert data into table fixed_pdb
'''
from helper import *
from bioomics import OpenmmPdbfixer

def retrieve_pdb(args, meta):
    pdb_group = args['pdb_group']
    table_name = args['table_name']
    query = f"""
        SELECT pdb_id, relative_path
        FROM meta_pdb
        WHERE pdb_group = {pdb_group}
            AND pdb_id NOT IN (
                SELECT pdb_id FROM {table_name}
            )
            AND pdb_id IN (
                SELECT pdb_id from view_antibody
            )
    ;"""
    qc = QueryComplex(args['verbose'])
    rows = qc.list_data(query)
    print(f"{len(rows)} PDB files ready...")
    meta['pdb_group'] = int(pdb_group)
    meta['rows'] += len(rows)
    return rows

def run_fixer(row, args, meta):
    time0 = datetime.now()
    pdb_id = row['pdb_id']
    pdb_endpoint = os.path.join(args['table_name'], pdb_id[:2])
    outdir = os.path.join(args['pdb_dir'], pdb_endpoint)
    output_prefix = os.path.join(outdir, pdb_id)
    
    fixed_pdb_file = output_prefix + '.pdb'
    if args['overwrite'] or not os.path.isfile(fixed_pdb_file):
        Dir(outdir).init_dir()
        # run pdbfixer
        input_pdb = os.path.join(args['pdb_dir'], row['relative_path'])
        op = OpenmmPdbfixer(input_pdb, args['verbose'])
        info = op.fast_fixer(output_prefix)
        meta[info['status']] += 1
    else:
        meta['skip'] += 1
    
    # build record
    rec = [pdb_id, None, time0]
    if os.path.isfile(fixed_pdb_file):
        rec[1] = os.path.join(pdb_endpoint, f"{pdb_id}.pdb")
    else:
        meta['null'] += 1
    return tuple(rec)


if __name__ == "__main__":
    args = {
        'verbose': True,
        'overwrite': False,
        'pdb_dir': os.getenv('pdb_dir'),
        'table_name': 'fixed_pdb',
    }
    args['pdb_group'] = sys.argv[1]

    rows = retrieve_pdb(args, meta)
    for row in rows:
        rec = run_fixer(row, args, meta)
        bc = BuildComplex([rec,], args['verbose'])
        bc.insert_batch_records(args['table_name'])

    print(meta)
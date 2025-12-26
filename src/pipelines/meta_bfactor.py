'''
requirements:
    - PDB raw data
    - table meta_pdb
example: 
    - parallel -j4 python src/pipelines/meta_bfactor.py ::: {0..116}
    - python src/pipelines/meta_bfactor.py 0
functions:
    - calculate B-factor of PDB
'''

from helper import *
from bioomics import ProcessPickle, CalculatePlddt

def retrieve_rows(pdb_group, args, meta):
    table_name = args['table_name']
    query = f"""
        SELECT pdb_id, relative_path
        FROM meta_pdb
        WHERE pdb_group = {pdb_group}
            AND pdb_id NOT IN (
                SELECT pdb_id FROM {table_name}
            );
    """
    # each record must have pdb_id and pdb_file
    rows = QueryComplex(args['verbose']).list_data(query)
    print(f"{len(rows)} PDB files ready...")
    meta['pdb_group'] = int(pdb_group)
    meta['rows'] += len(rows)
    for row in rows:
        pdb_file = os.path.join(args['pdb_dir'], row['relative_path'])
        if os.path.isfile(pdb_file):
            row['pdb_file'] =  pdb_file
            yield row

def cal_bfactor(rows, args, meta):
    for row in rows:
        pdb_id, pdb_file = row['pdb_id'], row['pdb_file']
        endpoint = os.path.join(args['table_name'], pdb_id[:2], pdb_id)
        outdir = os.path.join(args['data_dir'], endpoint)
        Dir(outdir).init_dir()

        # calculate bfactor
        file_name = f"bfactor.pkl"
        pfile = os.path.join(outdir, file_name)
        if args['overwrite'] or (not os.path.isfile(pfile)):
            # load structure before calculation
            cp = CalculatePlddt(pdb_file)
            cp.load_structure()
            
            if args['verbose']:
                print(f"Try to calclate B-factor given {pdb_file}...")
            data = cp.complex_ca_bfactor()
            if data.get('bfactor'):
                if args['verbose']:
                    print(f"Save B-factor results to {pfile}")
                status = ProcessPickle(pfile).save_pickle(data)
                if status:
                    meta['save_pickle'] += 1
                else:
                    meta['pickle_fail'] += 1
            else:
                meta['no_bfactor'] += 1
        else:
            meta['skip_bfactor'] += 1

        # build record
        rec = (pdb_id, None)
        if os.path.isfile(pfile):
            relative_path = os.path.join(endpoint, file_name)
            rec = (pdb_id, relative_path)
        yield rec


if __name__ == "__main__":
    args = {
        'verbose': True,
        'overwrite': False,
        'pdb_dir': os.getenv("pdb_dir"),
        'data_dir': os.getenv("data_dir"),
        'table_name': 'meta_bfactor',
    }
    pdb_group = sys.argv[1]
    row_iter = retrieve_rows(pdb_group, args, meta)

    # calculate B-factor
    record_iter = cal_bfactor(row_iter, args, meta)

    # put bfactor into table
    print(f"Try put meta data to {args['table_name']}...")
    bc = BuildComplex(record_iter, args['verbose'])
    bc.insert_batch_records(args['table_name'])

    print(meta)


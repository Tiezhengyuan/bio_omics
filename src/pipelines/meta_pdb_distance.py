'''
requirments:
    - PDB raw data
    - table meta_combo2, only calculate distance of two chains 
example: 
    - parallel -j8 python src/pipelines/meta_pdb_distance.py ::: {0..117}
    - python src/pipelines/meta_pdb_distance.py 0
functions:
    - calculate distance and analyze PPI
    - put meta data to table meta_pdb_distance
'''

from helper import *
from bioomics import DistanceCalculate, ProcessPickle

def retrieve_query(pdb_group, args, meta):
    table_name = args['table_name']
    query = f"""
        SELECT S.pdb_id, S.combo_id, S.chain_combo, S.combo_pdb
        FROM meta_combo2   S
        LEFT JOIN meta_pdb P ON S.pdb_id = P.pdb_id
        WHERE P.pdb_group = {pdb_group}
            AND P.relative_path IS NOT NULL
            AND S.combo_id IS NOT NULL
            AND S.combo_pdb IS NOT NULL
            AND S.combo_id NOT IN(
                SELECT combo_id FROM {table_name}
            );
    """
    rows = QueryComplex(args['verbose']).list_data(query)
    meta['rows'] = len(rows)
    print(f"Group {pdb_group}: Collect combo2 chains from {len(rows)} PDB ...")
    return rows

def cal_dist(rows, args, meta):
    for row in rows:
        pdb_id = row['pdb_id']
        endpoint = os.path.join(args['table_name'], pdb_id[:2], pdb_id)
        outdir = os.path.join(args['data_dir'], endpoint)
        Dir(outdir).init_dir()

        file_name =  f"{row['combo_id']}.pkl"
        pfile = os.path.join(outdir, file_name)
        if args['overwrite'] or (not os.path.isfile(pfile)):
            input_pdb = os.path.join(args['data_dir'], row['combo_pdb'])
            result = DistanceCalculate(input_pdb, args).process(row['chain_combo'])
            if result:
                result.update(row)
                ProcessPickle(pfile, args['verbose']).save_pickle(result)
                meta['dist_succeed'] += 1
            else:
                meta['dist_fail'] += 1
        else:
            meta['skip_cal_dist'] += 1
        
        # build record
        relative_path = os.path.join(endpoint, file_name) if os.path.isfile(pfile) else None
        yield pdb_id, row['combo_id'], relative_path

def validate_data(pdb_group, args, meta):
    table_name = args['table_name']
    uc = UpdateComplex(table_name, args['verbose'])

    query = f"""
        SELECT S.pdb_id, S.combo_id, S.relative_pkl
        FROM {table_name}   S
        LEFT JOIN meta_pdb P ON S.pdb_id = P.pdb_id
        WHERE P.pdb_group = {pdb_group}
            AND S.relative_pkl IS NOT NULL
    ;"""
    rows = QueryComplex(args['verbose']).list_data(query)
    for row in rows:
        pfile = os.path.join(args['data_dir'], row['relative_pkl'])
        res = ProcessPickle(pfile).is_good()
        if not res:
            msg = f"The file {pfile} is empty or wrong. {row}"
            uc.delete_rows('combo_id', row['combo_id'])
            msg += f" Delete {row} in {table_name}."
            meta['delete_rows'] += 1
            if os.path.isfile(pfile):
                msg += f" Try to delete {pfile}."
                os.remove(pfile)
                meta['delete_empty_pickle'] += 1
            print(msg)


if __name__ == "__main__":
    args = {
        'verbose': True,
        'overwrite': False,
        'chunk_size': 20,
        'data_dir': os.getenv("data_dir"),
        'table_name': 'meta_pdb_distance',
        'dist_cutoff': 10,
        'span': 10,
    }
    pdb_group = sys.argv[1]
    meta['pdb_group'] = pdb_group
    pdb_iter = retrieve_query(pdb_group, args, meta)
    # calculate atom distance
    record_iter = cal_dist(pdb_iter, args, meta)

    # put meta data to table
    bc = BuildComplex(record_iter, args['verbose'], args['chunk_size'])
    bc.insert_batch_records(args['table_name'])

    print("Try to validate pickle files...")
    validate_data(pdb_group, args, meta)

    print(meta)

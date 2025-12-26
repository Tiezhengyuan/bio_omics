'''
requirments:
    - PDB raw data
    - table meta_combo2_pdb
example: 
    - parallel -j4 python src/pipelines/pdb_sasa.py ::: {0..116}
    - python src/pipelines/pdb_sasa.py 0
functions:
    - calculate SASA and analyze PPI
    - put meta data to table meta_pdb_sasa
'''

from helper import *
from bioomics import SasaProcess, ProcessPickle

def retrieve_query(pdb_group, args, meta):
    table_name = args['table_name']
    # always return complex pdb
    query = f"""
        SELECT F.pdb_id, F.model_no, F.chain_combo, F.relative_path
        FROM meta_combo2_pdb F
        LEFT JOIN meta_pdb P ON F.pdb_id = P.pdb_id
        WHERE P.pdb_group = {pdb_group}
            AND F.relative_path IS NOT NULL
            AND F.pdb_id NOT IN(
                SELECT pdb_id FROM {table_name}
            );
    """
    df = QueryComplex(args['verbose']).list_data(query, True)
    if df is not None and len(df) > 0:
        meta['rows'] = len(df)
        g = df.groupby('pdb_id')
        return g
    return []

def get_pool(sub_rows):
    pool = []
    for _, row in sub_rows.iterrows():
        # bouned pdb file
        bfile = os.path.join(args['data_dir'], row['relative_path'])
        prefix = bfile[:bfile.rindex('-')-1]
        if os.path.isfile(bfile):
            unbounded_files = []
            chains = row['chain_combo'].split('-')
            for chain_no in chains:
                ufile = prefix + chain_no + '_0.pdb'
                if os.path.isfile(ufile):
                    unbounded_files.append(ufile)
            pool.append((bfile, unbounded_files))
    return pool

def cal_sasa(pdb_iter, args, meta):
    for pdb_id, sub in pdb_iter:
        endpoint = os.path.join(args['table_name'], pdb_id[:2], pdb_id)
        outdir = os.path.join(args['data_dir'], endpoint)
        Dir(outdir).init_dir()

        file_name =  'pair_sasa.pkl'
        pfile = os.path.join(outdir, file_name)
        if args['overwrite'] or (not os.path.isfile(pfile)):
            pdb_files_pool = get_pool(sub)
            cd = SasaProcess(verbose=args['verbose'])
            # result is list format, items by chain_combo
            data = cd.process(pdb_files_pool)
            if data:
                meta['pdb'] += 1
                meta['chain_combo'] += len(data)
                ProcessPickle(pfile).save_pickle(data)
                if args['verbose']:
                    print(f"Create a new pickle file {pfile}")
            else:
                meta['fail'] += 1
        else:
            meta['skip'] += 1
        
        # build record
        if os.path.isfile(pfile):
            relative_path = os.path.join(endpoint, file_name)
            yield pdb_id, relative_path
        else:
            yield pdb_id, None


if __name__ == "__main__":
    args = {
        'verbose': True,
        'overwrite': False,
        'chunk_size': 20,
        'data_dir': os.getenv('data_dir'),
        'table_name': 'meta_pdb_sasa',
    }
    pdb_group = sys.argv[1]
    pdb_iter = retrieve_query(pdb_group, args, meta)
    # calculate atom distance
    record_iter = cal_sasa(pdb_iter, args, meta)

    print(f"Try to put meta data to table {args['table_name']}")
    bc = BuildComplex(record_iter, args['verbose'], args['chunk_size'])
    bc.insert_batch_records(args['table_name'])

    print(meta)

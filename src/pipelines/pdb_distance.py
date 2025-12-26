'''
requirments:
    - PDB raw data
    - table meta_chain_pdb
example: 
    - parallel -j4 python src/pipelines/pdb_distance.py ::: {0..116}
    - python src/pipelines/pdb_distance.py 0
functions:
    - calculate distance and analyze PPI
    - put meta data to table meta_pdb_distance
'''

from helper import *
from bioomics import DistanceCalculate, ProcessPickle

def retrieve_query(pdb_group, args, meta):
    table_name = args['table_name']
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
    meta['rows'] = len(df)
    meta['pdb_group'] = int(pdb_group)
    if len(df) > 0:
        g = df.groupby('pdb_id')
        return g
    return []
    

def cal_dist(pdb_iter, args, meta):
    for pdb_id, sub in pdb_iter:
        endpoint = os.path.join(args['table_name'], pdb_id[:2], pdb_id)
        outdir = os.path.join(args['data_dir'], endpoint)
        Dir(outdir).init_dir()

        file_name = 'pair_distance.pkl'
        pfile = os.path.join(outdir, file_name)
        if args['overwrite'] or (not os.path.isfile(pfile)):
            data = []
            for _, row in sub.iterrows():
                input_pdb = os.path.join(args['data_dir'], row['relative_path'])
                cd = DistanceCalculate(input_pdb)
                item = cd.process(row['chain_combo'])
                if item:
                    item.update(row)
                    data.append(item)
                    meta['chain_combo'] += 1
                else:
                    meta['fail'] += 1
            if data:
                ProcessPickle(pfile).save_pickle(data)
        else:
            meta['skip'] += 1
        
        # build record
        if os.path.isfile(pfile):
            relative_path = os.path.join(endpoint, file_name)
            print(pdb_id, relative_path)
            yield pdb_id, relative_path
        else:
            yield pdb_id, None


if __name__ == "__main__":
    args = {
        'verbose': True,
        'overwrite': False,
        'chunk_size': 20,
        'data_dir': os.getenv("data_dir"),
        'table_name': 'meta_pdb_distance', 
    }
    pdb_group = sys.argv[1]
    pdb_iter = retrieve_query(pdb_group, args, meta)
     # calculate atom distance
    record_iter = cal_dist(pdb_iter, args, meta)

    # put meta data to table
    bc = BuildComplex(record_iter, args['verbose'], args['chunk_size'])
    bc.insert_batch_records(args['table_name'])

    print(meta)

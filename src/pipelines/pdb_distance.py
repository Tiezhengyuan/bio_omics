'''
requirments:
    - PDB raw data
    - table meta_chain_pdb
example: 
    - parallel -j4 python src/pipelines/pdb_distance.py ::: {0..46}
    - python src/pipelines/pdb_distance.py 46
functions:
    - calculate distance and analyze PPI
    - put meta data to table meta_pdb_distance
'''
from collections import defaultdict
from dotenv import load_dotenv
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import QueryComplex, BuildComplex, DistanceCalculate, ProcessPickle


def retrieve_query(pdb_group, args, meta):
    query = f"""
        SELECT F.pdb_id, F.model_no, F.chain_combo, F.pdb_file
        FROM meta_combo_pdb F
        LEFT JOIN meta_pdb P
        ON F.pdb_id = P.pdb_id
        WHERE P.pdb_group = {pdb_group}
            AND F.pdb_file IS NOT NULL
            AND F.pdb_id NOT IN(
                SELECT pdb_id
                FROM meta_pdb_distance
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
        pfile = os.path.join(args['data_dir'], 'data', pdb_id[:2], pdb_id, 'pair_distance.pkl')
        if args['overwrite'] or (not os.path.isfile(pfile)):
            data = []
            for _, row in sub.iterrows():
                cd = DistanceCalculate(row['pdb_file'])
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
        if os.path.isfile(pfile):
            yield pdb_id, pfile
        else:
            yield pdb_id, None


if __name__ == "__main__":
    load_dotenv()
    args = {
        'verbose': True,
        'overwrite': True,
        'data_dir': os.getenv("data_dir"),
        'chunk_size': 20,
    }
    meta = defaultdict(int)

    pdb_group = sys.argv[1]
    pdb_iter = retrieve_query(pdb_group, args, meta)
    # calculate atom distance
    record_iter = cal_dist(pdb_iter, args, meta)

    # put meta data to table meta_pdb_distance
    bc = BuildComplex(record_iter, args['verbose'], args['chunk_size'])
    bc.insert_meta_pdb_distance()

    print(meta)

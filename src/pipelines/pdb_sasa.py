'''
requirments:
    - PDB raw data
    - table meta_combo_pdb
example: 
    - parallel -j4 python src/pipelines/pdb_sasa.py ::: {0..46}
    - python src/pipelines/pdb_sasa.py 46
functions:
    - calculate SASA and analyze PPI
    - put meta data to table meta_pdb_sasa
'''
from collections import defaultdict
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import QueryComplex, BuildComplex, SasaProcess, ProcessPickle


def retrieve_query(pdb_group, args, meta):
    # always return complex pdb
    query = f"""
        SELECT F.pdb_id, F.model_no, F.chain_combo, F.pdb_file
        FROM meta_combo_pdb F
        LEFT JOIN meta_pdb P
        ON F.pdb_id = P.pdb_id
        WHERE P.pdb_group = {pdb_group}
            AND F.pdb_file IS NOT NULL
            AND F.pdb_id NOT IN(
                SELECT pdb_id
                FROM meta_pdb_sasa
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
        bfile = row['pdb_file']
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
        pfile = os.path.join(args['data_dir'], pdb_id[:2], pdb_id, 'pair_sasa.pkl')
        if args['overwrite'] or (not os.path.isfile(pfile)):
            pdb_files_pool = get_pool(sub)
            cd = SasaProcess(verbose=args['verbose'])
            # result is list format, items by chain_combo
            data = cd.process(pdb_files_pool)
            if data:
                meta['pdb'] += 1
                meta['chain_combo'] += len(data)
                ProcessPickle(pfile).save_pickle(data)
            else:
                meta['fail'] += 1
        else:
            meta['skip'] += 1
        # record iterator
        rec = (pdb_id, pfile) if os.path.isfile(pfile) else (pdb_id, None)
        yield rec

if __name__ == "__main__":
    pdb_group = sys.argv[1]
    args = {
        'verbose': True,
        'overwrite': True,
        'data_dir': '/home/yuan/data/pdb/data',
        'chunk_size': 20,
    }
    meta = defaultdict(int)
    pdb_iter = retrieve_query(pdb_group, args, meta)
    # calculate atom distance
    record_iter = cal_sasa(pdb_iter, args, meta)

    # put meta data to table meta_pdb_distance
    bc = BuildComplex(record_iter, args['verbose'], args['chunk_size'])
    bc.insert_meta_pdb_sasa()

    print(f"PDB group: {pdb_group}")
    print(meta)

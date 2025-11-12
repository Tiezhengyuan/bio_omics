'''
requirements:
    - PDB raw data
    - table pdb, meta_chain_pdb
example: 
    - parallel -j4 python src/pipelines/pdb_bfactor.py ::: {0..46}
    - python src/pipelines/pdb_bfactor.py 46
functions:
    - calculate B-factor by chains
    - put average bfactor to table pdb_bfactor
'''
from dotenv import load_dotenv
from collections import defaultdict
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import ProcessPickle, CalculatePlddt, QueryComplex, BuildComplex

def retrieve_rows(pdb_group, args, meta):
    query = f"""
        SELECT pdb_id, pdb_file
        FROM meta_pdb
        WHERE pdb_group = {pdb_group}
            AND pdb_id NOT IN (
                SELECT pdb_id
                FROM meta_bfactor
            );
    """
    # each record must have pdb_id and pdb_file
    rows = QueryComplex(args['verbose']).list_data(query)
    print(f"{len(rows)} PDB files ready...")
    meta['pdb_group'] = int(pdb_group)
    meta['rows'] += len(rows)
    return rows

def cal_bfactor(rows, args, meta):
    for row in rows:
        pdb_id = row['pdb_id']
        cascade_dir = os.path.join(args['outdir'], pdb_id[:2], pdb_id)
        pfile = os.path.join(cascade_dir, 'bfactor.pkl')
        if args['overwrite'] or (not os.path.isfile(pfile)):
            pdb_file = row['pdb_file']
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
        rec = (pdb_id, pfile) if os.path.isfile(pfile) else (pdb_id, None)
        yield rec


if __name__ == "__main__":
    meta = defaultdict(int)
    load_dotenv()
    args = {
        'verbose': True,
        'overwrite': False,
        'outdir': os.path.join(os.getenv("data_dir"), 'data'),
        'chunk_size': 20,
    }

    # retrieve pdb_ids
    pdb_group = sys.argv[1]
    rows = retrieve_rows(pdb_group, args, meta)

    # calculate B-factor
    record_iter = cal_bfactor(rows, args, meta)

    # put bfactor into table meta_bfactor
    print(f"Try put meta data to meta_bfactor...")
    bc = BuildComplex(record_iter, args['verbose'], args['chunk_size'])
    bc.insert_meta_bfactor()

    print(meta)


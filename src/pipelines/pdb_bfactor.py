'''
requirements:
    - PDB raw data
    - table pdb, meta_chain_pdb
example: 
    - parallel -j4 python /home/yuan/bio/bio_omics/src/pipelines/pdb_bfactor.py ::: {0..46}
    - python /home/yuan/bio/bio_omics/src/pipelines/pdb_bfactor.py 46
functions:
    - calculate B-factor by chains
    - put average bfactor to table pdb_bfactor
'''
from collections import defaultdict
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import ProcessPickle, CalculatePlddt, QueryComplex, BuildComplex

def retrieve_rows(pdb_group, meta):
    query = f"""
        SELECT pdb_id, pdb_file
        FROM meta_pdb
        WHERE pdb_group = {pdb_group}
            AND pdb_id NOT IN (
                SELECT DISTINCT pdb_id
                FROM pdb_bfactor
            );
    """
    # each record must have pdb_id and pdb_file
    rows = QueryComplex().list_data(query)
    print(f"{len(rows)} PDB files ready...")
    meta['rows'] += len(rows)
    return rows

def cal_bfactor(args, meta):
    pdb_id = args['pdb_id']
    cascade_dir = os.path.join(args['outdir'], pdb_id[:2], pdb_id)
    pfile = os.path.join(cascade_dir, 'bfactor.pkl')
    if args['overwrite'] or (not os.path.isfile(pfile)):
        # load structure before calculation
        cp = CalculatePlddt(args['pdb_file'], args['outdir'])
        cp.load_structure()
        if args['verbose']:
            print(f"Try to calclate B-factor given {args['pdb_file']}...")
        data = cp.complex_ca_bfactor()
        if data.get('bfactor'):
            if args['verbose']:
                print(f"Save B-factor results to {pfile}")
            status = ProcessPickle(pfile).save_pickle(data)
            meta[status] += 1
        else:
            meta['null'] += 1
    else:
        meta['skio'] += 1
    if os.path.isfile(pfile):
        return pfile

if __name__ == "__main__":
    meta = defaultdict(int)
    args = {
        'verbose': True,
        'overwrite': False,
        'outdir': '/home/yuan/data/pdb/data',
    }

    # retrieve pdb_ids
    pdb_group = sys.argv[1]
    rows = retrieve_rows(pdb_group, meta)

    # calculate B-factor
    pfiles = []
    for row in rows:
        row.update(args)
        pfile = cal_bfactor(row, meta)
        if pfile:
            rec = {
                'pdb_id': row['pdb_id'],
                'pfile': pfile,
            }
            pfiles.append(rec)

    # put bfactor into table pdb_bfactor
    if len(pfiles) > 0:
        print(f"Try retrieve data from {len(pfiles)}...")
        # insert data
        bc = BuildComplex(pfiles, args['verbose'])
        bc.insert_pdb_bfactor()
    else:
        print('No newly pickle files are detected. No insertion.')

    print(pdb_group, meta)


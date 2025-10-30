'''
requirements:
    - PDB raw data
    - table pdb, chain
example: 
    - parallel -j12 python /home/yuan/bio/bio_omics/src/pipelines/pdb_bfactor.py \
        ::: /home/yuan/data/pdb/_list/pdb_*.json
    - python /home/yuan/bio/bio_omics/src/pipelines/pdb_bfactor.py \
        /home/yuan/data/pdb/_list/pdb_46.json
functions:
    - calculate B-factor by chains
    - put average bfactor to table pdb_bfactor
'''
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import ProcessJson, ProcessPickle, CalculatePlddt, \
    QueryComplex, BuildComplex


def cal_bfactor(args):
    cp = CalculatePlddt(args['pdb_file'], args['outdir'])
    cascade_dir = cp.init_cascade_outdir(args['pdb_id'])
    pfile = os.path.join(cascade_dir, 'bfactor.pkl')
    if args['overwrite'] or (not os.path.isfile(pfile)):
        # load structure before calculation
        cp.load_structure()
        if args['verbose']:
            print(f"Try to calclate B-factor given {args['pdb_file']}...")
        data = cp.complex_ca_bfactor()
        if data.get('bfactor'):
            print(data)
            if args['verbose']:
                print(f"Save B-factor results to {pfile}")
            ProcessPickle(pfile).save_pickle(data)
            return 'succeed'
        else:
            return 'null'
    if os.path.isfile(pfile):
        args['pfile'] = pfile
    return 'skip'

def build_pdb_bfactor(args_pool):
    query = """
        SELECT pdb_id FROM pdb
        WHERE pdb_id NOT IN (
            SELECT DISTINCT pdb_id FROM pdb_bfactor
        );
    """
    rows = QueryComplex().list_data(query)
    un_pdb_ids = [i['pdb_id'] for i in rows]
    pfiles = []
    for args in args_pool:
        if 'pfile' in args and args['pdb_id'] in un_pdb_ids:
            pfiles.append(args)
    print(f"Try retrieve data from {len(pfiles)}...")

    # insert data
    bc = BuildComplex(pfiles, args['verbose'])
    bc.insert_pdb_bfactor()

if __name__ == "__main__":
    status = ['succeed', 'skip', 'null', 'na', 'fail']
    meta = {i:0 for i in status}
    json_file = sys.argv[1]
    args_pool = ProcessJson.load_json(json_file)

    # calculate B-factor
    for args in args_pool:
        args['verbose'] = True
        args['overwrite'] = False
        status = cal_bfactor(args)
        meta[status] += 1
    print(meta)

    # put bfactor into table pdb_bfactor
    build_pdb_bfactor(args_pool)
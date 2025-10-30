'''
requirements:
    - PDB raw data
    - table pdb, chain
example: 
    - parallel -j12 python /home/yuan/bio/bio_omics/src/pipelines/pdb_dssp.py \
        ::: /home/yuan/data/pdb/_list/pdb_*.json
    - python /home/yuan/bio/bio_omics/src/pipelines/pdb_dssp.py \
        /home/yuan/data/pdb/_list/pdb_46.json
functions:
    - Analyze secondary structure by chains
    - put dssp sequence to table pdb_dssp
'''
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import ProcessJson, QueryComplex, BuildComplex

def build_pdb_dssp(args, pdb_ids):
    # only consider model.id=0
    query = """
        SELECT pdb_id, chain_no FROM pdb_chain
        WHERE model_no = 0 AND pdb_id NOT IN (
            SELECT DISTINCT pdb_id FROM pdb_secondary
        );
    """
    rows = QueryComplex().list_data(query)
    files = []
    for row in rows:
        pdb_id, chain_no = row['pdb_id'], row['chain_no']
        if pdb_id in pdb_ids:
            pdb_file = os.path.join(args['data_dir'], pdb_id[:2], \
                pdb_id, 'chains', f'{pdb_id}_{chain_no}.pdb')
            if os.path.isfile(pdb_file):
                row['pdb_file'] = pdb_file
                files.append(row)
    print(f"Try retrieve data from {len(files)}...")

    # insert data
    bc = BuildComplex(files, args.get('verbose'))
    bc.insert_pdb_secondary()

if __name__ == "__main__":
    status = ['succeed', 'skip', 'null', 'na', 'fail']
    meta = {i:0 for i in status}
    args = {
        'verbose': True,
        'data_dir': '/home/yuan/data/pdb/data',
    }

    # load pdb_ids
    json_file = sys.argv[1]
    args_pool = ProcessJson.load_json(json_file)
    pdb_ids = {i['pdb_id']:i for i in args_pool}

    # put bfactor into table pdb_bfactor
    build_pdb_dssp(args, pdb_ids)
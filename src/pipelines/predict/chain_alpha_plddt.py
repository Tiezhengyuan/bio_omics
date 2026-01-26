'''
requirements:
    - alphafold prediciton model in PDB
    - table chain_alphafold
example: 
    - parallel -j8 python app.py chain_alpha_plddt ::: {0..117}
    - python app.py chain_alpha_plddt 0
functions:
    - calculate pLDDT
    - put pLDDT to table chain_alphafold_plddt
'''

from src.helper import *
from bioomics import CalculatePlddt

@header
def retrieve_rows(args, meta):
    pdb_group = args['pdb_group']
    table_name = args['table_name']
    query = f"""
        SELECT a.relative_pdb
        FROM {table_name}     a
        LEFT JOIN pdb_chainid c ON a.chain_id = c.chain_id
        LEFT JOIN meta_pdb    p ON c.pdb_id = p.pdb_id
        WHERE p.pdb_group = {pdb_group}
            AND a.relative_pdb IS NOT NULL
    ;"""
    curr_rows = QueryComplex().list_data(query)
    pdb_files = [i['relative_pdb'] for i in curr_rows]

    query = f"""
        SELECT c.pdb_id, a.chain_id, a.endpoint
        FROM chain_alphafold a
        LEFT JOIN pdb_chainid c ON a.chain_id = c.chain_id
        LEFT JOIN meta_pdb p ON c.pdb_id = p.pdb_id
        WHERE p.pdb_group = {pdb_group}
            AND a.endpoint IS NOT NULL
    ;"""
    rows = QueryComplex().list_data(query)
    meta['rows'] = len(rows)
    print(f"Filter {len(rows)} PDB data...")
    for row in rows:
        indir = os.path.join(args['predict_dir'], row['endpoint'])
        for file_name in os.listdir(indir):
            if file_name.endswith('.pdb'):
                relative_pdb = os.path.join(row['endpoint'], file_name)
                if relative_pdb not in pdb_files:
                    row['relative_pdb'] = relative_pdb
                    row['pdb_file'] = os.path.join(indir, file_name)
                    row['pkl_name'] = file_name.replace('.pdb', '') + '.pkl'
                    ranking = re.findall(r'rank_(\d+)_', file_name)
                    row['ranking'] = int(ranking[0]) if ranking else None
                    yield row

def cal_plddt(rows, args, meta):
    for row in rows:
        pdb_id = row['pdb_id']
        chain_id = row['chain_id']
        endpoint = os.path.join(args['table_name'], pdb_id[:2], pdb_id, chain_id)
        outdir = Path(args['predict_dir']) / endpoint
        outdir.mkdir(parents=True, exist_ok=True)

        data = None
        pfile = outdir / row['pkl_name']
        if args['overwrite'] or (not os.path.isfile(pfile)):
            # load structure before calculation
            cp = CalculatePlddt(row['pdb_file'], outdir)
            cp.load_structure()
            if args['verbose']:
                print(f"Try to calclate p given {row['pdb_file']}...")
            data = cp.complex_ca_plddt()
            if data.get('pLDDT'):
                if args['verbose']:
                    print(f"Save pLDDT results to {pfile}")
                status = ProcessPickle(pfile).save_pickle(data)
                meta[status] += 1
            else:
                meta['null'] += 1
        else:
            meta['skip'] += 1

        # build record
        if data and data.get('pLDDT'):
            relative_path = os.path.join(endpoint, row['pkl_name'])
            for item in data['pLDDT']:
                rec = (
                    row['chain_id'],
                    row['relative_pdb'],
                    row['ranking'],
                    'alphafold2',
                    item['stats'].get('mean_plddt'),
                    relative_path,
                )
                yield rec
 
if __name__ == "__main__":
    args = {
        'pdb_group': sys.argv[1],
        'verbose': True,
        'overwrite': False,
        'predict_dir': os.getenv('predict_dir'),
        'table_name': 'chain_plddt',
    }
    # retrieve pdb_ids
    rows = retrieve_rows(args, meta)

    # calculate pLDDT
    plddt_data_iter = cal_plddt(rows, args, meta)

    print("Try to insert pLDDT data into table")
    # insert data
    bc = BuildComplex(plddt_data_iter, args['verbose'])
    bc.insert_batch_records(args['table_name'])

    footer(meta)

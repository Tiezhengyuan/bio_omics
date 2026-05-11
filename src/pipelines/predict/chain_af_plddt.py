'''
requirements:
    - alphafold prediciton model in PDB
    - table chain_af
example: 
    - parallel -j8 python app.py chain_af_plddt ::: {0..240}
    - python app.py chain_af_plddt 0
functions:
    - calculate pLDDT
    - put pLDDT to table chain_af_plddt
'''
from src.helper import *
from bioomics import ProcessPickle, CalculatePlddt

@header
def retrieve_rows(params, meta):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
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
        FROM chain_af         a
        LEFT JOIN pdb_chainid c ON a.chain_id = c.chain_id
        LEFT JOIN meta_pdb    p ON c.pdb_id = p.pdb_id
        WHERE p.pdb_group = {pdb_group}
            AND a.endpoint IS NOT NULL
    ;"""
    rows = QueryComplex().list_data(query)
    meta['rows'] = len(rows)
    print(f"Begin to scan {len(rows)} PDB data...")
    for row in rows:
        indir = os.path.join(params['predict_dir'], row['endpoint'])
        for file_name in os.listdir(indir):
            if file_name.endswith('.pdb'):
                relative_pdb = os.path.join(row['endpoint'], file_name)
                # this file not in table chain_plddt
                if relative_pdb not in pdb_files:
                    row['relative_pdb'] = relative_pdb
                    row['pdb_file'] = os.path.join(indir, file_name)
                    row['pkl_name'] = file_name.replace('.pdb', '') + '.pkl'
                    ranking = re.findall(r'rank_(\d+)_', file_name)
                    row['ranking'] = int(ranking[0]) if ranking else None
                    yield row

def cal_plddt(rows, params, meta):
    for row in rows:
        meta['records'] += 1
        pdb_id = row['pdb_id']
        chain_id = row['chain_id']
        endpoint = os.path.join(params['table_name'], pdb_id[:2], pdb_id, chain_id)
        outdir = os.path.join(params['predict_dir'], endpoint)
        Path(outdir).mkdir(parents=True, exist_ok=True)

        pfile = os.path.join(outdir, row['pkl_name'])
        if params['overwrite'] or (not os.path.isfile(pfile)):
            # load structure before calculation
            cp = CalculatePlddt(row['pdb_file'], outdir)
            cp.load_structure()
            if params['verbose']:
                print(f"Try to calclate p given {row['pdb_file']}...")
            data = cp.complex_ca_plddt()
            if data.get('pLDDT'):
                if params['verbose']:
                    print(f"Save pLDDT results to {pfile}")
                status = ProcessPickle(pfile).save_pickle(data)
                meta[status] += 1
            else:
                meta['null'] += 1
        else:
            meta['skip'] += 1

        # build record
        data = ProcessPickle(pfile).load_pickle()
        if data and data.get('pLDDT'):
            meta['insert_chains'] += 1
            relative_path = os.path.join(endpoint, row['pkl_name'])
            for item in data['pLDDT']:
                meta['insert_records'] += 1
                rec = (
                    row['chain_id'],
                    row['relative_pdb'],
                    row['ranking'],
                    item['stats'].get('mean_plddt'),
                    relative_path,
                )
                yield rec
 
if __name__ == "__main__":
    params.update({
        'overwrite': False,
        'table_name': 'chain_af_plddt',
    })
    # retrieve pdb_ids
    rows = retrieve_rows(params, meta)

    # calculate pLDDT
    plddt_data_iter = cal_plddt(rows, params, meta)

    print("Try to insert pLDDT data into table")
    # insert data
    bc = BuildComplex(params['verbose'])
    bc.insert_batch_records(plddt_data_iter, params['table_name'])

    footer(meta)

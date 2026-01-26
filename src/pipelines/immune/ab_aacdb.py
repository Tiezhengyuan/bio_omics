'''
example: 
    - python src/pipelines/ab_aacdb.py
functions:
    - download data from AACDB
    - put data to table ab_aacdb
'''

from helper import *

def download(args):
    outdir = os.path.join(args['download_dir'], 'AACDB')
    Dir(outdir).init_dir()
    url = 'https://i.uestc.edu.cn/AACDB/data_zip'
    os.system(f"wget -c {url}/protein_table.txt -P {outdir}")

def retrieve_data(args, meta):
    infile = os.path.join(args['download_dir'], 'AACDB', 'protein_table.txt')
    df = pd.read_csv(infile, sep='\t')
    if args['verbose']:
        print('Size of AACDB: ', df.shape)
    meta['rows'] = len(df)
    rows = df.to_dict(orient='records')
    if args['verbose']:
        print('columns: ', list(rows[0]))
    for row in rows:
        pdb_id = row.get('pdb')
        if pdb_id:
            chains = row['chains'].split('_')
            record = (
                pdb_id.upper(),
                chains[0] if chains else None, #antiboby_chains
                chains[1] if len(chains) > 1 else None, # antigen_chains
                row.get('antibody'),
                row.get('protein'),
                row.get('INN(clinical_trial)'), #drug
                row.get('targets'),
                row.get('method'),
                row.get('organism'),
                row.get('resolution'),
                row.get('reference'),
            )
            yield record


if __name__ == "__main__":
    args = {
        'verbose': True,
        'chunk_size': 50,
        'download_dir': os.getenv("download_dir"),
        'table_name': 'ab_aacdb',
    }
    download(args)

    # prepare data
    record_iter = retrieve_data(args, meta)

    # insertion
    table_name = 'ab_aacdb'
    print(f"Try to update data from AACDB into table {table_name}.")
    bc = BuildComplex(record_iter, args['verbose'], args['chunk_size'])
    bc.empty_table(table_name)
    bc.insert_batch_records(args['table_name'])

    print(meta)




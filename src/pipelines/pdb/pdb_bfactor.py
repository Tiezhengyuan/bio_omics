'''
requirements:
    - bfactor.pkl
    - table meta_bfactor
example: 
    - parallel -j8 python app.py pdb_bfactor ::: {0..240}
    - python app.py pdb_bfactor 0
functions:
    - put average bfactor to table pdb_bfactor
'''

from src.helper import *
from bioomics import ProcessPickle

@header
def retrieve_pdb(params, meta):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        SELECT B.pdb_id, B.relative_pkl
        FROM meta_bfactor B, meta_pdb P
        WHERE B.pdb_id=P.pdb_id
            AND P.pdb_group = {pdb_group}
            AND B.relative_pkl IS NOT NULL
            AND B.pdb_id NOT IN (
                SELECT DISTINCT pdb_id FROM {table_name}
            );
    """
    # each record must have pdb_id and pdb_file
    rows = QueryComplex(params['verbose']).list_data(query)
    print(f"{len(rows)} files ready...")
    meta['rows'] += len(rows)
    for row in rows:
        row['pickle_file'] = os.path.join(params['output_pdb_dir'], row['relative_pkl'])
        if os.path.isfile(row['pickle_file']):
            yield row

def pdb_bfactor(row_iter, meta):
    for row in row_iter:
        tag = 0
        pdb_id = row['pdb_id']
        pfile = row['pickle_file']
        try:
            data = ProcessPickle(pfile).load_pickle()
            for item in data.get('bfactor', []):
                model_no = item['model'].id
                chain_no = item['chain'].id
                avg = item['stats'].get('mean')
                bfactor = item['bfactor']['bfactor'].astype(str).str.cat(sep=',')
                yield (pdb_id, model_no, chain_no, avg, bfactor)
                tag = 1
                meta['records'] += 1
        except Exception as e:
            meta['fail'] += 1
            if os.path.isfile(pfile):
                os.remove(pfile)
            print(f"Warning: pfile={pfile}, error={e}")
        if tag == 0:
            meta['null'] += 1
            yield (pdb_id,) + (None,) * 4


if __name__ == "__main__":
    params.update({
        'table_name': 'pdb_bfactor',
    })

    # retrieve B-factor
    row_iter = retrieve_pdb(params, meta)
    record_iter = pdb_bfactor(row_iter, meta)

    print(f"Try put meta data to {params['table_name']}...")
    bc = BuildComplex(params['verbose'])
    bc.insert_batch_records(record_iter, params['table_name'])

    footer(meta)
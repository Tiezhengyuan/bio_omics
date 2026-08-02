#!/usr/bin/python
'''
requirements:
    - PDB raw data should be downloaded to the local
    - table pdb
example:
    - bomics meta_pdb
function:
    - build table pdb and insert data into table meta_pdb
'''
from src.helper import *
from bioomics import ProcessJson

def get_meta_json(data, pdb_ids, params, meta):
    data_format = params['data_format']
    list_dir = params['output_pdb_dir'] / '_list' / data_format
    for file_name in os.listdir(list_dir):
        if file_name.endswith('.json'):
            pdb_group = os.path.splitext(file_name)[0]
            path = os.path.join(list_dir, file_name)
            print(f"Retrieve meta data from {path}")
            pdb_meta = ProcessJson.load_json(path)
            meta['pdb_group'] += 1
            for item in pdb_meta:
                pdb_id = item['pdb_id']
                if pdb_id not in pdb_ids:
                    relative_path = item['pdb_file'].replace(str(params['rawdata_dir']), '')
                    if relative_path.startswith('/'):
                        relative_path = relative_path[1:]
                    meta[f'new_{data_format}'] += 1
                    if pdb_id not in data:
                        data[pdb_id] = {
                            'pdb_id': pdb_id,
                            'pdb_group': int(pdb_group),
                        }
                    data[pdb_id][data_format] = relative_path
                else:
                    meta[f'skip_{data_format}'] += 1

def build(data, meta):
    for pdb_id, item in data.items():
        meta['insert_records'] += 1
        rec = (pdb_id, item.get('pdb'), item.get('mmCIF'), item['pdb_group'],)
        yield rec


def validate_meta_pdb(params, meta):
    table_name = params['table_name']
    print(f"Scan {table_name}...")
    query = f"""
        SELECT pdb_id, relative_pdb, relative_cif
        FROM {table_name}
    ;"""
    rows = QueryComplex(params['verbose']).list_data(query)
    del_pool = set()
    for row in rows:
        pdb_id = row['pdb_id']
        if row['relative_pdb']:
            pdb_file = params['rawdata_dir'] / row['relative_pdb']
            if not pdb_file.is_file():
                del_pool.add(pdb_id)
        if row['relative_cif']:
            pdb_file = params['rawdata_dir'] / row['relative_cif']
            if not pdb_file.is_file():
                del_pool.add(pdb_id)
    del_pool = list(del_pool)

    # delete rows
    n = len(del_pool)
    if n > 0:
        print(f"{n} out of {len(rows)} pdb/cif raw data are not existing.", end=' ')
        dc = DeleteComplex(params['verbose'])
        chunk_size = 100
        end = n - chunk_size if n > chunk_size else n
        for i in range(0, end, chunk_size):
            values = del_pool[i:i+chunk_size]
            res = dc.delete_rows(table_name, 'pdb_id', values)
            if res:
                meta['delete_rows'] += len(values)
                print(i, end=' ')
        print(f'{n}\n')
    else:
        print(f"PASS: Validate {table_name}.")


if __name__ == "__main__":
    params.update({
        'table_name': 'meta_pdb',
    })
    # validate existing pdb or cif files 
    validate_meta_pdb(params, meta)
    
    # retrieve pdb in database
    qc = QueryComplex(params['verbose'])
    pdb_ids = qc.unique_values(params['table_name'], 'pdb_id')

    data = {}
    for data_format in ('pdb', 'mmCIF'):
        params['data_format'] = data_format
        get_meta_json(data, pdb_ids, params, meta)
    record_iter = build(data, meta)

    print(f"Try to insert PDB meta data into table {params['table_name']}...")
    bc = BuildComplex(params['verbose'])
    bc.insert_batch_records(record_iter, params['table_name'])
    
    footer(params, meta)

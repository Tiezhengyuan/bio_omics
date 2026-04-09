'''
example:
    - parallel -j12 python app.py move ::: {0..237}
    - python app.py move 0
'''



from src.helper import *

def retrieve(params, meta):
    pdb_group = params['pdb_group']
    col_name = params['col_name']
    table_name = params['table_name']
    query = f"""
        SELECT a.pdb_id, a.{col_name}
        FROM {table_name}  a
        LEFT JOIN meta_pdb p ON a.pdb_id=p.pdb_id
        WHERE p.pdb_group = {pdb_group}
    ;"""
    rows = QueryComplex(params['verbose']).list_data(query)
    meta['rows'] = len(rows)

    source_dir = os.path.join(params['data_dir'], params['dir_name'])
    dest_dir = os.path.join(params['data_dir'], params['dir_name'] + '_new')
    for row in rows:
        pdb_id = row['pdb_id']
        dir1 = os.path.join(source_dir, pdb_id[:2], pdb_id, row[col_name])
        print(dir1)
        if os.path.exists(dir1):
            dir2 = os.path.join(dest_dir, pdb_id[:2], pdb_id, row[col_name])
            shutil.move(dir1, dir2)
            meta['move'] += 1
        else:
            meta['skip'] += 1

if __name__ == "__main__":
    params.update({
        'table_name': 'meta_combo2',
        'dir_name': 'combo2_relaxed',
        'col_name': 'combo_id',
    })
        
    retrieve(params, meta)
    footer(meta)



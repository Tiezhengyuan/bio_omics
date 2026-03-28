'''

example:
    - python app.py clean_predict 0
    - parallel -j8 python app.py clean_predict ::: {0..117}
'''
from src.helper import *
import shutil


def not_existing_endpoint(params, table_name):
    '''
    endpoint is db but no local directory then delete the db rows
    '''
    pdb_group = params['pdb_group']
    query = f"""
        select a.pdb_id, a.endpoint
        from {table_name} a
        left join meta_pdb p on a.pdb_id = p.pdb_id
        where a.endpoint is not null
            and p.pdb_group = {pdb_group}
    ;"""
    rows = QueryComplex(True).list_data(query)
    meta[f'{table_name}_endpoints'] = len(rows)

    ud = UpdateComplex(table_name, True)
    for row in rows:
        _dir = os.path.join(params['predict_dir'], row['endpoint'])
        # delete rows if directory is not existings
        if not os.path.isdir(_dir):
            if not params['dry_run']:
                print(f"Try to delete record {row['endpoint']} from {table_name}")
                ud.delete_rows('endpoint', [row['endpoint'],])
                meta[f'{table_name}_delete_records'] += 1
            else:
                print(f"Recommend deletion record: {row['endpoint']}")
                meta[f'{table_name}_recommend_delete_record'] += 1

def delete_uncomplete_prediction(params, table_name, rows):
    '''
    there is local dir but 
    '''
    meta[f'{table_name}_rows'] = len(rows)
    for row in rows:
        pdb_id, _id = row['pdb_id'], row['id']
        outdir = os.path.join(params['predict_dir'], table_name, pdb_id[:2], pdb_id, _id)
        if os.path.isdir(outdir):
            meta[f'{table_name}_outdir'] += 1
            done_txt = os.path.join(outdir, f"{_id}.done.txt")
            # delete directory if no dont_text
            if not os.path.isfile(done_txt):
                if not params['dry_run']:
                    print(f"Try to delete directory {outdir}")
                    shutil.rmtree(outdir)
                    meta[f'{table_name}_delete_dir'] += 1
                else:
                    print(f"Recommend deletion directory: {outdir}")
                    meta[f'{table_name}_recommend_delete_dir'] += 1

def chain_prediction(params, table_name):
    pdb_group = params['pdb_group']
    query = f"""
        select c.pdb_id, c.chain_id as id
        from pdb_chainid c
        left join meta_pdb p on c.pdb_id = p.pdb_id
        where p.pdb_group = {pdb_group}
    ;"""
    rows = QueryComplex(True).list_data(query)
    delete_uncomplete_prediction(params, table_name, rows)


def combo2_prediction(params, table_name):
    pdb_group = params['pdb_group']
    query = f"""
        select c.pdb_id, c.combo_id as id
        from meta_combo2 c
        left join meta_pdb p on c.pdb_id = p.pdb_id
        where p.pdb_group = {pdb_group}
    ;"""
    rows = QueryComplex(True).list_data(query)
    delete_uncomplete_prediction(params, table_name, rows)




if __name__ == "__main__":
    params = {
        'pdb_group': sys.argv[1],
        'data_dir': os.getenv('data_dir'),
        'predict_dir': os.getenv('predict_dir'),
        # show attempts but not actual deletion
        'dry_run': True,
    }

    table_name = 'chain_alphafold'
    # delete uncomplete prediction
    # chain_prediction(params, table_name)
    # remove record of which endpoint is not existing
    # not_existing_endpoint(params, table_name)

    table_name =  'combo2_alphafold'
    # delete uncomplete prediction
    combo2_prediction(params, table_name)
    # remove record of which endpoint is not existing
    not_existing_endpoint(params, table_name)

    footer(meta)

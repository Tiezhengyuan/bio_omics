'''
example:
    python app.py collect_ab_distance
'''
from src.helper import *

def get_meta(table_name):
    query = f"""
        select D.pdb_id, D.combo_id, C.chain_combo, D.relative_pkl
        from {table_name} D
        left join meta_combo2 C on D.combo_id = C.combo_id
        where D.pdb_id in (
            select distinct pdb_id from view_antibody
        )
        AND D.relative_pkl is not null
    ;"""
    rows = QueryComplex(True).list_data(query)
    df = pd.DataFrame(rows)
    info = {
        'table_name': table_name,
        'pdb': df['pdb_id'].nunique(),
        'combo': df['combo_id'].nunique(),
        'regions': len(rows),
    }
    print(info)
    return rows

def abag_distance(rows, abag, params, meta):
    '''
    distance of antibody in antibody-antigen dimer
    '''
    print('Number of distance', len(rows))
    table_name = params['table_name']
    data = []
    for row in rows:
        combo_id = row['combo_id']
        if combo_id in abag:
            ab_chain_no = abag[combo_id]['ab_chain_no']
            rec = {
                'pdb_id': row['pdb_id'],
                'combo_id': combo_id,
                'chain_combo': row['chain_combo'],
                'ab_chain_no': ab_chain_no,
            }
            pfile = os.path.join(params['data_dir'], row['relative_pkl'])
            try:
                data = ProcessPickle(pfile).load_pickle()
                if data['chain1'] == ab_chain_no:
                    rec['dist'] = data['dist1']
                    data.append(rec)
                    meta[f'{table_name}_dist1'] += 1
                elif data['chain2'] == ab_chain_no:
                    rec['dist'] = data['dist2']
                    data.append(rec)
                    meta[f'{table_name}_dist2'] += 1
            except Exception as e:
                meta[f'{table_name}_fail'] += 1
                print(e)
    return data
    
if __name__ == "__main__":
    params = {
        'verbose': True,
        'overwrite': True,
        'data_dir': os.getenv('data_dir'),
    }
    outdir = Path(params['data_dir']) / 'collect' / 'distance'
    outdir.mkdir(parents=True, exist_ok=True)
    
    # antibody-antigen dimer in dict
    abag = QueryComplex(params['verbose']).abag_combo2()

    tables = ['combo2_distance', 'pdb_distance']
    for table_name in tables:
        rows = get_meta(table_name)
        # distance of antibody in antibody-antigen dimer
        outfile = outdir / f"distance_abag_{table_name}.pkl"
        data = {}
        if params['overwrite'] or not outfile.is_file():
            print(f"Try to count contacts of Ca: {table_name}...")
            params['table_name'] = table_name
            data = abag_distance(rows, abag, params, meta)
            if data is not None:
                print(f"Save {outfile}")
                ProcessPickle(outfile).save_pickle(data)


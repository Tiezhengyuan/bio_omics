'''
requirments:
    - table predict_summary
example: 
    - parallel -j8 python app.py predict_sasa ::: {0..10}
    - python app.py predict_sasa 0
functions:
    - calculate SASA and analyze PPI
    - put meta data to table predict_sasa
'''
from src.helper import *
from bioomics import SasaProcess, ProcessPickle, TABLE_COLUMNS

@header
def retrieve(params, meta):
    min_id = int(params['pdb_group']) * 1_000
    max_id = min_id + 1_000
    design_name = params['design_name']
    table_name = params['table_name']
    col_name = params['col_name']
    query = f"""
        SELECT a.predict_id, a.summary_id, a.{col_name} AS relative_pdb, d.pdb_files
        FROM predict_summary    a
        LEFT JOIN predict       b ON a.predict_id = b.predict_id
        LEFT JOIN {design_name} d ON b.design_id = d.design_id
        WHERE a.predict_id >= {min_id} AND a.predict_id < {max_id}
            AND a.{col_name} IS NOT NULL
            AND a.summary_id NOT IN (
                SELECT summary_id FROM {table_name}
            )
    """
    rows = QueryComplex(params['verbose']).list_data(query)
    meta['rows'] = len(rows)
    for row in rows:
        pdb_file = os.path.join(params['design_dir'], row['relative_pdb'])
        if os.path.isfile(pdb_file):
            row['model_pdb_file'] = pdb_file
            yield row

def cal_sasa(row_iter, params, meta):
    for row in row_iter:
        predict_id = str(row['predict_id'])
        endpoint = os.path.join(params['table_name'], predict_id[:2], predict_id)
        output_dir = os.path.join(params['design_dir'], endpoint)
        Path(output_dir).mkdir(parents=True, exist_ok=True)

        file_name = re.sub('pdb$', 'pkl', os.path.basename(row['relative_pdb']))
        pfile = os.path.join(output_dir, file_name)
        if params['overwrite'] or (not os.path.isfile(pfile)):
            try:
                pdb_files_pool = [
                    (row['model_pdb_file'], row['pdb_files'].split(',')),
                ]
                print(pdb_files_pool)
                SasaProcess(verbose=params['verbose']).process(pdb_files_pool, pfile)
                meta['sasa_succeed'] += 1
            except Exception as e:
                print(e)
                meta['sasa_fail'] += 1
        else:
            meta['skip'] += 1
        
        # build record
        if os.path.isfile(pfile):
            row['relative_pkl'] = os.path.join(endpoint, file_name)
            data = ProcessPickle(pfile, params['verbose']).load_pickle([])
            if data:
                row['delta_total_sasa'] = data[0].get('delta_total_sasa')
        rec = [row.get(col) for col in TABLE_COLUMNS[params['table_name']]]
        print(rec)
        yield tuple(rec)


if __name__ == "__main__":
    params.update({
        'overwrite': False,
        'chunk_size': 20,
    })
    for design_name in ('design_dimer', 'design_trimer'):
        params['design_name'] = design_name
        for name in ('unrelaxed', 'relaxed'):
            params['table_name'] = f'{name}_sasa'
            params['col_name'] = f"{name}_pdb"
            pdb_iter = retrieve(params, meta)
            record_iter = cal_sasa(pdb_iter, params, meta)

            bc = BuildComplex(params['verbose'], params['chunk_size'])
            bc.insert_batch_records(record_iter, params['table_name'])

    footer(meta)

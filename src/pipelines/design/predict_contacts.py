#!/usr/bin/python
'''
requirements:
    - predict, predict_summary
example:
    - parallel -j8 python app.py predict_contacts ::: {0..10}
    - python app.py predict_contacts 0
function:
    - simulate binding dissociation rate
    - insert data into table predict_contacts
'''
from src.helper import *
from bioomics import EvalPred, TABLE_COLUMNS

@header
def retrieve(params, meta):
    min_id = int(params['pdb_group']) * 1_000
    max_id = min_id + 1_000
    table_name = params['table_name']
    col_name = params['col_name']
    query = f"""
        SELECT summary_id, predict_id, {col_name} AS relative_pdb
        FROM predict_summary
        WHERE predict_id >= {min_id} AND predict_id < {max_id}
            AND {col_name} IS NOT NULL
            AND summary_id NOT IN (
                SELECT summary_id from {table_name}
            )
    """
    rows = QueryComplex().list_data(query)
    meta['rows'] = len(rows)
    for row in rows:
        pdb_file = os.path.join(params['design_dir'], row['relative_pdb'])
        if os.path.isfile(pdb_file):
            row['pdb_file'] = pdb_file
            yield row

def cal_contacts(rows, params, meta):
    cols = TABLE_COLUMNS[params['table_name']]
    for row in rows:
        predict_id = str(row['predict_id'])
        endpoint = os.path.join(params['table_name'], predict_id[:2], predict_id)
        output_dir = os.path.join(params['design_dir'], endpoint)
        Path(output_dir).mkdir(parents=True, exist_ok=True)

        # run prodigy
        file_name = re.sub('pdb$', 'prodigy', os.path.basename(row['relative_pdb']))
        outfile = os.path.join(output_dir, file_name)
        if params['overwrite'] or not os.path.isfile(outfile):
            EvalPred().prodigy(row['pdb_file'], outfile)
            meta['run_prodigy'] += 1

        # retrieve data
        if os.path.isfile(outfile):
            # retrieve data
            with open(outfile, 'r') as f:
                data = f.read()
            meta['output'] += 1
            result = EvalPred().parse_prodigy(data)
            row['binding_affinity'] = result.get('binding_affinity')
            row['dissociation_constant'] = result.get('dissociation_constant')
            row['relative_prodigy'] = os.path.join(endpoint, file_name)
        else:
            meta['no_output'] += 1

        # build record
        rec = [row.get(col) for col in cols]
        if params['verbose']:
            print(rec)
        yield tuple(rec)

if __name__ == "__main__":
    params.update({
        'overwrite': False,
        'chunk_size': 50,
    })

    for name in ('unrelaxed', 'relaxed'):
        params['table_name'] = f'{name}_contacts'
        params['col_name'] = f"{name}_pdb"
        rows = retrieve(params, meta)
        record_iter = cal_contacts(rows, params, meta,)
        bc = BuildComplex(params['verbose'], params['chunk_size'])
        bc.insert_batch_records(record_iter, params['table_name'])

    footer(meta)
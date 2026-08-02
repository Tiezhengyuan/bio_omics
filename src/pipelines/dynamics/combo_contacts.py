#!/usr/bin/python
'''
requirements:
    - table combo2, combo3, combo4, combo5
example:
    - parallel -j8 bomics combo_contacts ::: {0..241}
    - bomics combo_contacts 0
function:
    - simulate binding dissociation rate
    - insert data into table combo_contacts
'''
from src.helper import *
from bioomics import EvalPred, TABLE_COLUMNS

@header
def retrieve_pdb(params, meta):
    pdb_group = params['pdb_group']
    combo_name = params['combo_name']
    table_name = params['table_name']
    query = f"""
        SELECT DISTINCT c.pdb_id, c.combo_id, c.combo_pdb
        FROM label_complex     l
        LEFT JOIN {combo_name} c ON l.combo_id = c.combo_id
        LEFT JOIN meta_pdb     p ON c.pdb_id = p.pdb_id
        WHERE p.pdb_group = {pdb_group}
            AND NOT EXISTS (
                SELECT 1 FROM {table_name} b
                WHERE b.combo_id = c.combo_id
            )
    """
    rows = QueryComplex(params['verbose']).list_data(query)
    print(f"{len(rows)} PDB files ready...")
    meta[f'{combo_name}_rows'] += len(rows)
    return rows

def cal_contacts(rows, params, meta):
    combo_name = params['combo_name']
    for row in rows:
        pdb_id, combo_id = row['pdb_id'], row['combo_id']
        endpoint = Path(params['table_name']) / pdb_id[:2] / pdb_id
        output_dir = params['simulate_dir'] / endpoint
        output_dir.mkdir(parents=True, exist_ok=True)

        # run prodigy
        file_name = f'{combo_id}.txt'
        outfile = output_dir / file_name
        if params['overwrite'] or not os.path.isfile(outfile):
            input_pdb = params['output_pdb_dir'] / row['combo_pdb']
            EvalPred().prodigy(input_pdb, outfile)
            meta[f'{combo_name}_run_prodigy'] += 1

        # build record
        item = {
            'pdb_id': pdb_id,
            'combo_id': combo_id,
        }
        if os.path.isfile(outfile):
            # retrieve data
            with open(outfile, 'r') as f:
                data = f.read()
            meta[f'{combo_name}_output'] += 1
            result = EvalPred().parse_prodigy(data)
            item['binding_affinity'] = result.get('binding_affinity')
            item['dissociation_constant'] = result.get('dissociation_constant')
            item['relative_prodigy'] = str(endpoint / file_name)
        else:
            meta[f'{combo_name}_no_output'] += 1
        rec = [item.get(col) for col in TABLE_COLUMNS[params['table_name']]]
        if params['verbose']:
            print(rec)
        yield tuple(rec)


if __name__ == "__main__":
    params.update({
        'overwrite': False,
        'chunk_size': 50,
        'table_name': 'combo_contacts',
    })
    # clean null data in the table
    validate_data('relative_prodigy', params, meta, params['simulate_dir'])
    
    # cover all combo complex
    for combo_name in ('combo2', 'combo3', 'combo4', 'combo5'):
        params['combo_name'] = combo_name
        rows = retrieve_pdb(params, meta)
        record_iter = cal_contacts(rows, params, meta,)
        bc = BuildComplex(params['verbose'], params['chunk_size'])
        bc.insert_batch_records(record_iter, params['table_name'])

    # validate_data('relative_prodigy', params, meta)
    footer(meta)
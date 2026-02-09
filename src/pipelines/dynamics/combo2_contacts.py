#!/usr/bin/python
'''
requirements:
    - table meta_combo2
example:
    - parallel -j8 python app.py combo2_contacts ::: {0..117}
    - python app.py combo2_contacts 0
function:
    - fix pdb files
    - relax pdb
    - insert data into table combo2_contact
'''
from src.helper import *
from bioomics import EvalPred

@header
def retrieve_pdb(args, meta):
    pdb_group = args['pdb_group']
    table_name = args['table_name']
    query = f"""
        SELECT C.pdb_id, C.combo_id, C.combo_pdb
        FROM meta_combo2   C
        LEFT JOIN meta_pdb P ON C.pdb_id = P.pdb_id
        WHERE P.pdb_group = {pdb_group}
            AND C.combo_id NOT IN (
                SELECT combo_id FROM {table_name}
            )
    ;"""
    rows = QueryComplex(args['verbose']).list_data(query)
    print(f"{len(rows)} PDB files ready...")
    meta['rows'] += len(rows)
    return rows

def cal_contacts(rows, args, meta):
    for row in rows:
        pdb_id, combo_id = row['pdb_id'], row['combo_id']
        endpoint = os.path.join(args['table_name'], pdb_id[:2], pdb_id)
        output_dir = os.path.join(args['data_dir'], endpoint)
        Dir(output_dir).init_dir()

        # run prodigy
        file_name = f'{combo_id}.txt'
        outfile = os.path.join(output_dir, file_name)
        if args['overwrite'] or not os.path.isfile(outfile):
            input_pdb = os.path.join(args['data_dir'], row['combo_pdb'])
            EvalPred().prodigy(input_pdb, outfile)
            meta['run_prodigy'] += 1

        # build record
        if os.path.isfile(outfile):
            # retrieve data
            with open(outfile, 'r') as f:
                data = f.read()
            meta['has_output'] += 1
            result = EvalPred().parse_prodigy(data)
            aff = result.get('binding_affinity')
            kd = result.get('dissociation_constant')
            path = os.path.join(endpoint, file_name)
            yield (pdb_id, combo_id, aff, kd, path)
        else:
            meta['no_output'] += 1
            yield (pdb_id, combo_id, None, None, None)


if __name__ == "__main__":
    args = {
        'pdb_group': sys.argv[1],
        'verbose': True,
        'overwrite': False,
        'chunk_size': 50,
        'data_dir': os.getenv('data_dir'),
        'table_name': 'combo2_contacts',
    }
    rows = retrieve_pdb(args, meta)
    record_iter = cal_contacts(rows, args, meta,)
    bc = BuildComplex(record_iter, args['verbose'], args['chunk_size'])
    bc.insert_batch_records(args['table_name'])

    validate_data('relative_prodigy', args, meta)
    footer(meta)
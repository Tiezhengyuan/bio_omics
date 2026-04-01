#!/usr/bin/python
'''
requirements:
    - table meta_combo3
example:
    - parallel --delay 5 -j8 python app.py combo3_relaxed ::: {0..99}
    - python app.py combo3_relaxed
function:
    - fix pdb files
    - relax pdb
    - insert data into table combo3_relaxed
'''
from src.helper import *
from bioomics import OpenmmPdbfixer, OpenmmSimulation, ProcessJson

@header
def retrieve_pdb(params, meta):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        SELECT c.pdb_id, c.combo_id, c.combo_pdb
        FROM meta_combo3   c
        LEFT JOIN meta_pdb p ON c.pdb_id = c.pdb_id
        WHERE p.pdb_group = {pdb_group}
            AND c.combo_id NOT IN (
                SELECT combo_id FROM {table_name}
            )
            AND c.pdb_id IN (
                SELECT pdb_id from view_antibody
            )
    ;"""
    rows = QueryComplex(params['verbose']).list_data(query)
    print(f"{len(rows)} PDB files ready...")
    meta['rows'] += len(rows)
    return rows

def fix_and_relax(row, params, meta):
    time0 = datetime.now()
    pdb_id, combo_id = row['pdb_id'], row['combo_id']
    endpoint = os.path.join(params['table_name'], pdb_id[:2], pdb_id, combo_id)
    output_dir = os.path.join(params['data_dir'], endpoint)

    with WriteTmp(params['tmp_dir'], output_dir, params['table_name']) as writer:
        outdir = writer.temp_dir_path
        info = {
            'pdb_id': pdb_id,
            'combo_id': combo_id,
            'endpoint': endpoint,
            'output_dir': output_dir,
            'outdir': outdir,
            'time0': time0.strftime("%Y-%m-%d %H:%M:%S"),
        }

        # fix pdb
        fixed_pdb = os.path.join(outdir, 'fixed.pdb')
        if params['overwrite'] or not os.path.isfile(fixed_pdb):
            input_pdb = os.path.join(params['data_dir'], row['combo_pdb'])
            curr_params = {
                'verbose': params['verbose'],
                'outdir': outdir,
            }
            op = OpenmmPdbfixer(input_pdb, curr_params)
            fix_info = op.fast_fixer()
            info['input_pdb_file'] = input_pdb
            info['fix'] = fix_info
            meta['fix_' + fix_info['status']] += 1
        else:
            fix_info_json = os.path.join(outdir, 'fixed_info.json')
            info['fix'] = ProcessJson.load_json(fix_info_json)
            meta['fix_skip'] += 1

        # relax pdb
        if os.path.isfile(fixed_pdb):
            info['fixed_pdb_file'] = fixed_pdb
            info['fixed_pdb'] = os.path.join(endpoint, "fixed.pdb")
            # simulate pdb
            relaxed_pdb = os.path.join(outdir, 'relaxed.pdb')
            if params['overwrite'] or not os.path.isfile(relaxed_pdb):
                curr_params = {
                    'pdb_id': pdb_id,
                    'outdir': outdir,
                }
                curr_params.update(params)
                relax_info = OpenmmSimulation(fixed_pdb, curr_params).relax_pdb()
                info['relax'] = relax_info
                meta['relax_' + relax_info['status']] += 1
            else:
                relax_info_json = os.path.join(outdir, 'relaxed_info.json')
                info['relax'] = ProcessJson.load_json(relax_info_json)
                meta['relax_skip'] += 1                    
            # update record
            if os.path.isfile(relaxed_pdb):
                info['relaxed_pdb_file'] = relaxed_pdb
                info['relaxed_pdb'] = os.path.join(endpoint, "relaxed.pdb")
        else:
            meta['relax_null'] += 1
        
        # save info
        info_json = os.path.join(outdir, 'info.json')
        if not os.path.isfile(info_json):
            info['time1'] = datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            ProcessJson.save_json(info, info_json)

    # update record
    rec = (pdb_id, combo_id, time0,
        info.get('fixed_pdb'), info.get('relaxed_pdb'),
        info['relax'].get('original_energy') if 'relax' in info else None,
        info['relax'].get('relaxed_energy') if 'relax' in info else None,
    )
    if params['verbose']:
        print('Record: ', rec)
    return rec


if __name__ == "__main__":
    params = {
        'verbose': True,
        'overwrite': False,
        'data_dir': os.getenv('data_dir'),
        'tmp_dir': os.getenv('tmp_dir'),
        'table_name': 'combo3_relaxed',
        'force_fields': ['amber14-all.xml', 'amber14/tip3pfb.xml'],
        'max_iterations': 2000,
        'steps': 10_000,
    }
    add_params(params)

    rows = retrieve_pdb(params, meta)
    for row in rows:
        rec = fix_and_relax(row, params, meta,)
        bc = BuildComplex([rec,], params['verbose'])
        bc.insert_batch_records(params['table_name'])

    validate_data('fixed_pdb', params, meta)
    validate_data('relaxed_pdb', params, meta)
    footer(params, meta)
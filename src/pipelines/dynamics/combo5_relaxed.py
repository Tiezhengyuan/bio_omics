#!/usr/bin/python
'''
requirements:
    - table meta_combo4
example:
    - parallel -j12 python app.py combo5_relaxed ::: {0..240}
    - python app.py combo5_relaxed 0
function:
    - fix pdb files
    - relax pdb
    - insert data into table combo5_relaxed
'''
from src.helper import *
from bioomics import WriteTmp, ProcessJson, OpenmmPdbfixer, OpenmmSimulation

@header
def retrieve_pdb(params, meta):
    q = QueryComplex(params['verbose'])
    rows = q.retrieve_combo(params['pdb_group'], params['combo'], params['table_name'])
    print(f"{len(rows)} PDB files ready...")
    meta['rows'] += len(rows)
    return rows

def dry_run(rows, params, meta):
    for row in rows:
        time0 = datetime.now()
        pdb_id, combo_id = row['pdb_id'], row['combo_id']
        endpoint = os.path.join(params['table_name'], pdb_id[:2], pdb_id, combo_id)
        output_dir = os.path.join(params['simulate_dir'], endpoint)

        print("Try to rebuild info.json...")
        info = {
            'pdb_id': pdb_id,
            'combo_id': combo_id,
            'endpoint': endpoint,
            'output_dir': output_dir,
            'time0': time0.strftime("%Y-%m-%d %H:%M:%S"),
            'input_pdb_file': os.path.join(params['output_pdb_dir'], row['combo_pdb'])
        }
        # fix meta
        fixed_pdb = os.path.join(output_dir, 'fixed.pdb')
        if os.path.isfile(fixed_pdb):
            info['fixed_pdb_file'] = fixed_pdb
            info['fixed_pdb'] = os.path.join(endpoint, "fixed.pdb")
            fix_info_json = os.path.join(output_dir, 'fixed_info.json')
            info['fix'] = ProcessJson.load_json(fix_info_json)
        # relax meta
        relaxed_pdb = os.path.join(output_dir, 'relaxed.pdb')
        if os.path.isfile(relaxed_pdb):
            info['relaxed_pdb_file'] = relaxed_pdb
            info['relaxed_pdb'] = os.path.join(endpoint, "relaxed.pdb")
            relax_info_json = os.path.join(output_dir, 'relaxed_info.json')
            info['relax'] = ProcessJson.load_json(relax_info_json)
        info_json = os.path.join(output_dir, 'info.json')
        ProcessJson.save_json(info, info_json)

        # build record for insertion if data is ready
        if info and 'fixed_pdb' in info:
            # update record
            rec = (
                pdb_id, combo_id, info['time0'], info['fixed_pdb'], info.get('relaxed_pdb'),
                info['relax'].get('original_energy') if 'relax' in info else None,
                info['relax'].get('relaxed_energy') if 'relax' in info else None,
            )
            meta['dry_run'] += 1
            if params['verbose']:
                print('Record: ', rec)
            yield rec

def fix_and_relax(row, params, meta):
    time0 = datetime.now()
    pdb_id, combo_id = row['pdb_id'], row['combo_id']
    endpoint = os.path.join(params['table_name'], pdb_id[:2], pdb_id, combo_id)
    output_dir = os.path.join(params['simulate_dir'], endpoint)

    with WriteTmp(params['tmp_dir'], output_dir, params['table_name']) as writer:
        outdir = writer.temp_dir_path
        info = {
            'pdb_id': pdb_id,
            'combo_id': combo_id,
            'endpoint': endpoint,
            'output_dir': output_dir,
            'time0': time0.strftime("%Y-%m-%d %H:%M:%S"),
        }

        # fix pdb
        fixed_pdb = os.path.join(outdir, 'fixed.pdb')
        if params['overwrite'] is True or not os.path.isfile(fixed_pdb):
            input_pdb = os.path.join(params['output_pdb_dir'], row['combo_pdb'])
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
            meta['fix_skip'] += 1
        # load fix_info
        fix_info_json = os.path.join(outdir, 'fixed_info.json')
        if os.path.isfile(fix_info_json):
            info['fix'] = ProcessJson.load_json(fix_info_json)

        # relax pdb
        relaxed_pdb = os.path.join(outdir, 'relaxed.pdb')
        if os.path.isfile(fixed_pdb):
            info['fixed_pdb_file'] = fixed_pdb
            info['fixed_pdb'] = os.path.join(endpoint, "fixed.pdb")
            # simulate pdb
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
                meta['relax_skip'] += 1                    
            # update record
            if os.path.isfile(relaxed_pdb):
                info['relaxed_pdb_file'] = relaxed_pdb
                info['relaxed_pdb'] = os.path.join(endpoint, "relaxed.pdb")
        else:
            meta['relax_null'] += 1
        # load relax info
        relax_info_json = os.path.join(outdir, 'relaxed_info.json')
        if os.path.isfile(relax_info_json):
            info['relax'] = ProcessJson.load_json(relax_info_json)        
        
        # save info
        info_json = os.path.join(outdir, 'info.json')
        info['time1'] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        ProcessJson.save_json(info, info_json)
    
    # upload results if needed
    copy_data(output_dir)

    # update record
    if 'fixed_pdb' in info:
        rec = (
            pdb_id, combo_id, time0, info['fixed_pdb'], info.get('relaxed_pdb'),
            info['relax'].get('original_energy') if 'relax' in info else None,
            info['relax'].get('relaxed_energy') if 'relax' in info else None,
        )
        if params['verbose']:
            print('Record: ', rec)
        meta['insertion'] += 1
        return rec

if __name__ == "__main__":
    params.update({
        'dry_run': False,
        'overwrite': False,
        'tmp_dir': os.getenv('tmp_dir'),
        'combo': 5,
        'table_name': 'combo5_relaxed',
        'force_fields': ['amber14-all.xml', 'amber14/tip3pfb.xml'],
        'max_iterations': 2_000,
        'steps': 10_000,
    })
    rows = retrieve_pdb(params, meta)

    if params['dry_run']:
        record_iter = dry_run(rows, params, meta)
        bc = BuildComplex(params['verbose'], 2)
        bc.insert_batch_records(record_iter, params['table_name'])
    else:
        for row in rows:
            rec = fix_and_relax(row, params, meta)
            bc = BuildComplex(params['verbose'])
            bc.insert_batch_records([rec,], params['table_name'])

    # validate_data('fixed_pdb', params, meta)
    # validate_data('relaxed_pdb', params, meta)
    footer(meta)
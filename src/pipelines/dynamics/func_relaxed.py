#!/usr/bin/python
'''
define funcitons used by the below
    - combo2_relaxed.py
    - combo3_relaxed.py
    - combo4_relaxed.py
    - combo5_relaxed.py
'''
from src.helper import *
from bioomics import WriteTmp, ProcessJson, OpenmmPdbfixer, OpenmmSimulation

@header
def retrieve_pdb(params, meta):
    # clean null data in the table
    validate_data('relaxed_pdb', params, meta, params['simulate_dir'])

    q = QueryComplex(params['verbose'])
    rows = q.retrieve_combo(params['pdb_group'], params['combo'], params['table_name'])
    meta['rows'] += len(rows)
    return rows

def dry_run(rows, params, meta):
    for row in rows:
        print(row)
        time0 = datetime.now()
        pdb_id, combo_id = row['pdb_id'], row['combo_id']
        endpoint = Path(params['table_name']) / pdb_id[:2] / pdb_id / combo_id
        output_dir = params['simulate_dir'] / endpoint

        print("Try to rebuild info.json...")
        info = {
            'pdb_id': pdb_id,
            'combo_id': combo_id,
            'endpoint': str(endpoint),
            'output_dir': str(output_dir),
            'time0': time0.strftime("%Y-%m-%d %H:%M:%S"),
            'input_pdb_file': str(params['output_pdb_dir'] / row['combo_pdb'])
        }
        # fix meta
        fixed_pdb = output_dir / 'fixed.pdb'
        if fixed_pdb.is_file():
            info['fixed_pdb_file'] = str(fixed_pdb)
            info['fixed_pdb'] = str(endpoint / "fixed.pdb")
            fix_info_json = output_dir / 'fixed_info.json'
            info['fix'] = ProcessJson.load_json(fix_info_json)
        # relax meta
        relaxed_pdb = output_dir / 'relaxed.pdb'
        if relaxed_pdb.is_file():
            info['relaxed_pdb_file'] = str(relaxed_pdb)
            info['relaxed_pdb'] = str(endpoint / "relaxed.pdb")
            relax_info_json = output_dir / 'relaxed_info.json'
            info['relax'] = ProcessJson.load_json(relax_info_json)
        info_json = output_dir / 'info.json'
        ProcessJson.save_json(info, info_json)

        # build record for insertion if data is ready
        if info and 'fixed_pdb' in info:
            # update record
            rec = (
                pdb_id,
                combo_id,
                info['time0'],
                info['fixed_pdb'],
                info.get('relaxed_pdb'),
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
    endpoint = Path(params['table_name']) / pdb_id[:2] / pdb_id / combo_id
    output_dir = params['simulate_dir'] / endpoint

    with WriteTmp(params['tmp_dir'], output_dir, params['table_name']) as writer:
        outdir = writer.temp_dir_path
        info = {
            'pdb_id': pdb_id,
            'combo_id': combo_id,
            'endpoint': str(endpoint),
            'output_dir': str(output_dir),
            'time0': time0.strftime("%Y-%m-%d %H:%M:%S"),
        }

        # fix pdb
        fixed_pdb = outdir / 'fixed.pdb'
        if params['overwrite'] is True or not fixed_pdb.is_file():
            input_pdb = params['output_pdb_dir'] / row['combo_pdb']
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
        fix_info_json = outdir / 'fixed_info.json'
        if fix_info_json.is_file():
            info['fix'] = ProcessJson.load_json(fix_info_json)


        # relax pdb
        relaxed_pdb = outdir / 'relaxed.pdb'
        if fixed_pdb.is_file():
            info['fixed_pdb_file'] = fixed_pdb
            info['fixed_pdb'] = str(endpoint / "fixed.pdb")
            # simulate pdb
            if params['overwrite'] or not relaxed_pdb.is_file():
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
            if relaxed_pdb.is_file():
                info['relaxed_pdb_file'] = str(relaxed_pdb)
                info['relaxed_pdb'] = str(endpoint / "relaxed.pdb")
        else:
            meta['relax_null'] += 1
        # load relax info
        relax_info_json = outdir / 'relaxed_info.json'
        if relax_info_json.is_file():
            info['relax'] = ProcessJson.load_json(relax_info_json)        
        
        # save info
        info_json = outdir / 'info.json'
        info['time1'] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        ProcessJson.save_json(info, info_json)
    
    # upload results if needed
    # depreciated in the future: both pdb and simulate directory are nfs sharing in remote host
    # copy_data(output_dir)

    # update record
    if info and 'fixed_pdb' in info:
        rec = (
            pdb_id,
            combo_id,
            time0,
            info['fixed_pdb'],
            info.get('relaxed_pdb'),
            info['relax'].get('original_energy') if 'relax' in info else None,
            info['relax'].get('relaxed_energy') if 'relax' in info else None,
        )
        if params['verbose']:
            print('Record: ', rec)
        meta['insertion'] += 1
        return rec
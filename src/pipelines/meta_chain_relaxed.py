#!/usr/bin/python
'''
requirements:
    - table meta_chain_pdb
example:
    - parallel -j4 python src/pipelines/meta_chain_relaxed.py ::: {0..116}
    - python src/pipelines/meta_chain_relaxed.py 0
function:
    - fix pdb files
    - relax pdb
    - insert data into table meta_chain_relaxed
'''
from helper import *
from bioomics import OpenmmPdbfixer, OpenmmSimulation, ProcessJson

def retrieve_pdb(args, meta):
    pdb_group = args['pdb_group']
    table_name = args['table_name']
    query = f"""
        SELECT C.pdb_id, C.chain_id, C.relative_path
        FROM meta_chain_pdb C
        LEFT JOIN pdb_chainid   I ON C.chain_id = I.chain_id
        LEFT JOIN pdb_chain_seq S ON C.chain_id = S.chain_id
        LEFT JOIN meta_pdb P ON C.pdb_id = P.pdb_id
        WHERE P.pdb_group = {pdb_group}
            AND I.model_no = 0
            AND S.chain_seq IS NOT NULL
            AND C.chain_id NOT IN (
                SELECT chain_id FROM {table_name}
            )
            AND C.pdb_id IN (
                SELECT pdb_id from view_antibody
            )
    ;"""
    qc = QueryComplex(args['verbose'])
    rows = qc.list_data(query)
    print(f"{len(rows)} PDB files ready...")
    meta['pdb_group'] = int(pdb_group)
    meta['rows'] += len(rows)
    return rows

def fix_and_relax(row, args, meta):
    time0 = datetime.now()
    pdb_id, chain_id = row['pdb_id'], row['chain_id']
    endpoint = os.path.join(args['table_name'], pdb_id[:2], pdb_id, chain_id)
    outdir = os.path.join(args['data_dir'], endpoint)
    info = {
        'pdb_id': pdb_id,
        'chain_id': chain_id,
        'endpoint': endpoint,
        'outdir': outdir,
        'time0': time0.strftime("%Y-%m-%d %H:%M:%S"),
    }

    # fix pdb
    fixed_pdb = os.path.join(outdir, 'fixed.pdb')
    if args['overwrite'] or not os.path.isfile(fixed_pdb):
        input_pdb = os.path.join(args['data_dir'], row['relative_path'])
        curr_args = {
            'verbose': args['verbose'],
            'outdir': outdir,
        }
        op = OpenmmPdbfixer(input_pdb, curr_args)
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
        if args['overwrite'] or not os.path.isfile(relaxed_pdb):
            curr_args = {
                'pdb_id': pdb_id,
                'outdir': outdir,
            }
            curr_args.update(args)
            relax_info = OpenmmSimulation(fixed_pdb, curr_args).relax_pdb()
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
    rec = (pdb_id, chain_id, time0,
        info.get('fixed_pdb'), info.get('relaxed_pdb'),
        info['relax'].get('original_energy') if 'relax' in info else None,
        info['relax'].get('relaxed_energy') if 'relax' in info else None,
    )
    if args['verbose']:
        print('Record: ', rec)
    return rec


if __name__ == "__main__":
    args = {
        'verbose': True,
        'overwrite': False,
        'data_dir': os.getenv('data_dir'),
        'table_name': 'meta_chain_relaxed',
        'force_fields': ['amber14-all.xml', 'amber14/tip3pfb.xml'],
        'max_iterations': 5000,
        'steps': 100_000,
    }
    args['pdb_group'] = sys.argv[1]

    rows = retrieve_pdb(args, meta)
    for row in rows:
        rec = fix_and_relax(row, args, meta)
        bc = BuildComplex([rec,], args['verbose'])
        bc.insert_batch_records(args['table_name'])

    print(meta)
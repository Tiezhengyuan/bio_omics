'''
requirements:
    - table af_summary
example: 
    - parallel -j8 bomics af_plddt ::: {0..241}
    - bomics af_plddt 0
functions:
    - calculate pLDDT by chain given prediction pdb
'''
from src.helper import *
from bioomics import ProcessPickle, CalculatePlddt


@header
def retrieve_rows(params, meta):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        SELECT a.pdb_id, a.category, a.input_id, a.ranking, a.relaxed_pdb
        FROM af_summary    a
        LEFT JOIN meta_pdb p ON a.pdb_id = p.pdb_id
        WHERE p.pdb_group = {pdb_group}
            AND a.relaxed_pdb IS NOT NULL
            AND a.relaxed_pdb NOT IN (
                SELECT relative_pdb FROM {table_name}
            )
    ;"""
    rows = QueryComplex().list_data(query)
    meta['raws'] = len(rows)
    for row in rows:
        pdb_file = params['predict_dir'] / row['relaxed_pdb']
        if pdb_file.is_file():
            row['pdb_file'] = pdb_file
            row['pkl_file_name'] = pdb_file.stem + '.pkl'
            meta['relaxed_pdb'] = len(rows)
            yield row

def cal_plddt(rows, params, meta):
    for row in rows:
        pdb_id = row['pdb_id']
        endpoint = Path(params['table_name']) / pdb_id[:2] / pdb_id / row['input_id']
        outdir = params['predict_dir'] / endpoint
        outdir.mkdir(parents=True, exist_ok=True)

        pfile = outdir / row['pkl_file_name']
        if params['overwrite'] or not pfile.is_file():
            # load structure before calculation
            cp = CalculatePlddt(row['pdb_file'], outdir)
            cp.load_structure()
            if params['verbose']:
                print(f"Try to calclate pLDDT by chain with the input {row['pdb_file']}...")
            data = cp.complex_ca_plddt()
            if data.get('pLDDT'):
                if params['verbose']:
                    print(f"Save pLDDT results to {pfile}")
                status = ProcessPickle(pfile).save_pickle(data)
                meta[status] += 1
            else:
                meta['null'] += 1
        else:
            meta['skip'] += 1

        # build record
        data = ProcessPickle(pfile).load_pickle()
        if data and data.get('pLDDT'):
            meta['insert_chains'] += 1
            for item in data['pLDDT']:
                meta['insert_records'] += 1
                rec = (
                    row['pdb_id'],
                    row['input_id'],
                    row['relaxed_pdb'],
                    str(endpoint / row['pkl_file_name']),
                    item['model'].id,
                    item['chain'].id,
                    item['stats'].get('mean_plddt'),
                )
                yield rec

 
if __name__ == "__main__":
    params.update({
        'overwrite': False,
        'table_name': 'af_plddt',
    })
        
    # retrieve pdb_ids
    rows = retrieve_rows(params, meta)

    # calculate pLDDT
    plddt_data_iter = cal_plddt(rows, params, meta)

    # insert data
    bc = BuildComplex(params['verbose'])
    bc.insert_batch_records(plddt_data_iter, params['table_name'])

    footer(meta)

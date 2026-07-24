'''
requirements:
    - alphafold prediciton model in PDB
    - table chain_af
example: 
    - parallel -j8 python app.py predict_plddt ::: {0..10}
    - python app.py predict_plddt 0
functions:
    - calculate pLDDT
    - put pLDDT to table predict_plddt
'''
from src.helper import *
from bioomics import ProcessPickle, CalculatePlddt

@header
def retrieve(params, meta):
    min_id = int(params['pdb_group']) * 1_000
    max_id = min_id + 1_000
    table_name = params['table_name']
    col_name = params['col_name']

    query = f"""
        SELECT predict_id, summary_id, {col_name} AS relative_pdb
        FROM predict_summary
        WHERE predict_id >= {min_id} AND predict_id < {max_id}
            AND {col_name} IS NOT NULL
            AND summary_id NOT IN (
                SELECT summary_id FROM {table_name}
            )
    """
    rows = QueryComplex().list_data(query)
    meta['rows'] = len(rows)
    for row in rows:
        pdb_file = os.path.join(params['design_dir'], row['relative_pdb'])
        if os.path.isfile(pdb_file):
            relative_pkl = re.sub('.pdb$', '.pkl', row['relative_pdb'])
            row['pdb_file'] = pdb_file,
            row['relative_pkl'] = relative_pkl
            yield row

def cal_plddt(rows, params, meta):
    for row in rows:
        meta['records'] += 1
        pfile = os.path.join(params['design_dir'], row['relative_pkl'])
        if params['overwrite'] or (not os.path.isfile(pfile)):
            # load structure before calculation
            cp = CalculatePlddt(row['pdb_file'])
            cp.load_structure()
            if params['verbose']:
                print(f"Try to calclate pLDDT given {row['pdb_file']}...")
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
                    row['summary_id'],
                    row['relative_pdb'],
                    row['relative_pkl'],
                    item['chain'].id,
                    item['stats'].get('mean_plddt'),
                )
                yield rec
 
if __name__ == "__main__":
    params.update({
        'overwrite': False,
        'chunk_size': 50,
    })
    for name in ('unrelaxed', 'relaxed'):
        params['table_name'] = f'{name}_plddt'
        params['col_name'] = f"{name}_pdb"
        # retrieve pdb_ids
        rows = retrieve(params, meta)
        record_iter = cal_plddt(rows, params, meta)

        print("Try to insert pLDDT data into table")
        bc = BuildComplex(params['verbose'], params['chunk_size'])
        bc.insert_batch_records(record_iter, params['table_name'])

    footer(meta)

'''
requirements:
    - PDB raw data
    - table meta_pdb
example: 
    - parallel -j8 python app.py meta_bfactor ::: {0..240}
    - python app.py meta_bfactor 0
functions:
    - calculate B-factor of PDB
'''

from src.helper import *
from bioomics import ProcessPickle, CalculatePlddt

@header
def retrieve_rows(params, meta):
    qc = QueryComplex(params['verbose'])
    rows = qc.meta_pdb(params['pdb_group'], params['table_name'])
    print(f"{len(rows)} PDB files ready...")
    meta['rows'] += len(rows)
    for row in rows:
        pdb_file = os.path.join(params['rawdata_pdb_dir'], row['relative_pdb'])
        if os.path.isfile(pdb_file):
            row['pdb_file'] =  pdb_file
            yield row

def cal_bfactor(rows, params, meta):
    for row in rows:
        pdb_id, pdb_file = row['pdb_id'], row['pdb_file']
        endpoint = os.path.join(params['table_name'], pdb_id[:2], pdb_id)
        outdir = os.path.join(params['output_pdb_dir'], endpoint)
        Path(outdir).mkdir(parents=True, exist_ok=True)

        # calculate bfactor
        file_name = f"bfactor.pkl"
        pfile = os.path.join(outdir, file_name)
        if params['overwrite'] or (not os.path.isfile(pfile)):
            # load structure before calculation
            cp = CalculatePlddt(pdb_file)
            cp.load_structure()
            
            if params['verbose']:
                print(f"Try to calclate B-factor given {pdb_file}...")
            data = cp.complex_ca_bfactor()
            if data.get('bfactor'):
                if params['verbose']:
                    print(f"Save B-factor results to {pfile}")
                status = ProcessPickle(pfile).save_pickle(data)
                if status:
                    meta['save_pickle'] += 1
                else:
                    meta['pickle_fail'] += 1
            else:
                meta['no_bfactor'] += 1
        else:
            meta['skip_bfactor'] += 1

        # build record
        rec = (pdb_id, None)
        if os.path.isfile(pfile):
            relative_path = os.path.join(endpoint, file_name)
            rec = (pdb_id, relative_path)
        yield rec


if __name__ == "__main__":
    params.update({
        'overwrite': False,
        'table_name': 'meta_bfactor',
    })
    row_iter = retrieve_rows(params, meta)

    # calculate B-factor
    record_iter = cal_bfactor(row_iter, params, meta)

    # put bfactor into table
    print(f"Try put meta data to {params['table_name']}...")
    bc = BuildComplex(params['verbose'])
    bc.insert_batch_records(record_iter, params['table_name'])

    footer(meta)

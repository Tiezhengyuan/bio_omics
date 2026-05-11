'''
requirements:
    - PDB files should be synchronized to the local.
    - The table pdb is ready.
example:
    - parallel -j8 python app.py pdb_compound ::: {0..240}
    - python app.py pdb_compound 0
function:
    - retrieve compound in head of PDB
'''
from src.helper import *
from bioomics import Pdb

def retrieve(params, meta):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        SELECT pdb_id, relative_path
        FROM meta_pdb
        WHERE pdb_group = {pdb_group}
            AND pdb_id NOT IN (
                SELECT DISTINCT pdb_id FROM {table_name}
            );
    """
    rows = QueryComplex(params['verbose']).list_data(query)
    meta['pdb_group'] = int(pdb_group)
    meta['rows'] += len(rows)
    for row in rows:
        row['pdb_file'] = os.path.join(params['rawdata_pdb_dir'], row['relative_path'])
        yield row

def build(rows, params, meta):
    for row in rows:
        pdb_id = row['pdb_id']
        pdb_file = row['pdb_file']
        # load pdb
        p = Pdb(pdb_file, None, params['verbose'])
        p.load_structure()
        # build record
        for compound_no, cp in p.compound.items():
            compound_no = int(compound_no)
            _text = []
            for _key, _value in cp.items():
                if _key not in ('chain',) and _value:
                    _text.append(f"{_key}: {_value}")
            compound_text = ' ; '.join(_text)
            for chain_no in cp.get('chain').split(','):
                chain_no = chain_no.strip()
                rec = (pdb_id, compound_no, chain_no, compound_text)
                meta['records'] += 1
                yield rec

##################################
if __name__ == "__main__":
    params.update({
        'chunk_size': 200,
        'table_name': 'pdb_compound',
    })

    # retrieve data
    rows = retrieve(params, meta)
    record_iter = build(rows, params, meta)

    print(f"Try to insert PDB data into table {params['table_name']}...")
    bc = BuildComplex(params['verbose'], params['chunk_size'])
    bc.insert_batch_records(record_iter, params['table_name'])

    footer(meta)

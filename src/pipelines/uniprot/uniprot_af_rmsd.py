'''
requirements:
    - PDB: meta_pdb, pdb_chain
    - blastp of pdb chain and table align_uniprot
    - alphafold prediction in PDB
example: 
    - parallel -j8 python app.py uniprot_af_rmsd ::: {0..240}
    - python app.py uniprot_af_rmsd 0
functions:
    - similarity between model and sample structure
'''

from src.helper import *
from bioomics import EvalPred

@header
def search_pdb(params, meta):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        SELECT C.chain_id, C.pdb_id, U.uniprot_acc, C.chain_pdb, A.pdb_path
        FROM chain_pdb          C
        LEFT JOIN meta_pdb      M ON C.pdb_id = M.pdb_id
        LEFT JOIN pdb_chain_seq S ON C.chain_id = S.chain_id
        LEFT JOIN uniprot       U ON S.chain_seq = U.pro_seq
        LEFT JOIN uniprot_af    A ON U.uniprot_acc = A.uniprot_acc
        WHERE M.pdb_group = {pdb_group}
            AND C.chain_pdb IS NOT NULL
            AND A.pdb_path IS NOT NULL
            AND C.chain_id NOT IN (
                SELECT chain_id FROM {table_name}
            );
    """
    rows = QueryComplex(params['verbose']).list_data(query)
    meta['rows'] = len(rows)
    for row in rows:
        row['native_pdb'] = os.path.join(params['output_pdb_dir'], row['chain_pdb'])
        row['model_pdb'] = os.path.join(params['rawdata_pdb_dir'], row['pdb_path'])
        if os.path.isfile(row['native_pdb']) and os.path.isfile(row['model_pdb']):
            yield row

def calculate_rmsd(rows, meta):
    for row in rows:
        rmsd = None
        try:
            rmsd = EvalPred().rmsd(row['native_pdb'], row['model_pdb'])
            if rmsd:
                meta['succeed_rmsd'] += 1
            else:
                meta['null_rmsd'] += 1    
        except Exception as e:
            meta['fail_rmsd'] += 1
        # wrap record
        rec = (row['chain_id'], row['uniprot_acc'], rmsd,\
            row['chain_pdb'], row['pdb_path'])
        yield rec


if __name__ == "__main__":
    params.update({
        'chunk_size': 10,
        'table_name': 'uniprot_af_rmsd',
    })

    # retrieve input
    row_iter = search_pdb(params, meta)
    # RMSD
    record_iter = calculate_rmsd(row_iter, meta)

    # insert data
    bc = BuildComplex(params['verbose'], params['chunk_size'])
    bc.insert_batch_records(record_iter, params['table_name'])

    footer(meta)

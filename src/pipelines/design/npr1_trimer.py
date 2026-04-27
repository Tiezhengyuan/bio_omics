#!/usr/bin/python
'''
requirements:
    - table view_antibody, label_chain, lablel_complex
example:
    - python app.py npr1_trimer
function:
    - design complex for prediction
'''
from src.helper import *
from bioomics import ProcessSeq, ProcessFasta

def get_npr1() -> tuple:
    # 9BCN is NPR-1 dimer of A-B chains
    chain_id = '9BCN_0_A_0'
    query = f"""
        select chain_seq, chain_pdb from view_pdb_chain
        where chain_id = '{chain_id}'
    """
    rows = QueryComplex(True).list_data(query)
    meta['ab_chain'] = len(rows)
    return chain_id, rows[0]['chain_seq'], rows[0]['chain_pdb']

# antibody chain
def get_complex(params):
    category = 'HLAg'
    # antigen
    target_id, target_seq, target_pdb = get_npr1()
    # antibody
    query = """
        SELECT a.pdb_id, a.combo_id, a.chain_combo, a.chain_label,
            a.chain_ids, b.combo_pdb, b.combo_fasta
        FROM label_complex a, combo2 b
        WHERE a.combo_id = b.combo_id AND a.category = 'AbAb'
    ;"""
    rows = QueryComplex(True).list_data(query)
    
    # remove duplicated design_id
    res = {} 
    for row in rows:
        row['target_id'] = target_id
        row['binder_id'] = row['combo_id']
        design_id = f"{target_id}-{row['binder_id']}"
        combo_fasta = os.path.join(params['data_dir'], row['combo_fasta'])
        binder_seq = ProcessFasta(combo_fasta).get_first_seq()
        if binder_seq:
            row['combo_seqs'] = [target_seq,] + binder_seq.split(':')
            row['chain_labels'] = f"{params['target_name']}-{row['chain_label']}"
            row['binder_pdb'] = os.path.join(params['data_dir'], row['combo_pdb'])
            row['target_pdb'] = os.path.join(params['data_dir'], target_pdb)
            row['endpoint'] = os.path.join(params['table_name'], row['chain_labels'])
            row['category'] = category
            res[design_id] = row
    return category, res

#create possible dimers
def design(rows:dict, params, meta):
    for design_id, row in rows.items():
        # prepare outdir
        outdir = os.path.join(params['data_dir'], row['endpoint'])
        Path(outdir).mkdir(parents=True, exist_ok=True)
        # prepare query_fasta
        fa_file_name = f"{design_id}.faa"
        fa_file = os.path.join(outdir, fa_file_name)
        if params['overwrite'] or not os.path.isfile(fa_file):
            combo_seq = {design_id: row['combo_seqs'],}
            status = ProcessSeq.to_colabfold_fasta(combo_seq, fa_file)
            meta['combo_fasta_'+ status] += 1
        else:
            meta['combo_fasta_skip'] += 1

        # prepare record
        if os.path.exists(fa_file):
            combo = 3
            yield (
                design_id, row['category'], row['chain_labels'], combo,
                params['target_name'], row['target_id'], row['binder_id'],
                fa_file, row.get('target_pdb'), row.get('binder_pdb'),
            )
        else:
            print(f"Skip {design_id}: {fa_file} is not existing")

if __name__ == "__main__":
    params.update({
        'overwrite': True,
        'chunk_size': 50,
        'table_name': 'design_multimer',
        'target_name': 'NPR1',
    })
    # retrieve data
    category, rows = get_complex(params)
    # delete rows before insertion
    UpdateComplex(params['table_name'], True).delete_rows('category', [category,])
    # build records
    record_iter = design(rows, params, meta)
    # insert records
    bc = BuildComplex(record_iter, params['verbose'], params['chunk_size'])
    bc.insert_batch_records(params['table_name'])

    footer(params, meta)
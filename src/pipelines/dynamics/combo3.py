'''
requirements:
    - The table meta_pdb is ready.
    - meta_chain_pdb is ready
example:
    - parallel -j8 python app.py combo3 ::: {0..240}
    - python app.py combo3 0
function:
    - split PDB of complex to paired chains PDB.
    - export complex csv and fasta for structure predication
note:
    - only consider model.id == 0
'''
from src.helper import *
from bioomics import Pdb, ProcessSeq


@header
def retrieve_data(params, meta):
    pdb_group = params['pdb_group']

    c = QueryComplex(params['verbose'])
    #combo_ids
    chains_num = params['combo']
    table_name = params['table_name']
    query = f"""
        SELECT c.pdb_id, c.model_no, c.combo_id, c.chain_combo,
            c.chain_label, c.chain_ids, p.relative_path
        FROM label_complex   c
        LEFT JOIN meta_pdb p on c.pdb_id = p.pdb_id
        WHERE p.pdb_group = {pdb_group}
            AND combo = {chains_num}
            AND c.combo_id NOT IN (
                SELECT combo_id from {table_name}
            )
    ;"""
    rows = c.list_data(query)
    ids = [i['combo_id'] for i in rows]
    print('Check duplicates:', len(rows), len(set(ids)))
    meta[f'combo_{chains_num}'] = len(rows)

    # chain_ids
    pdb_chains = c.pdb_chains(pdb_group)
    # add key: input_pdb, chain_seqs, chain_ids
    for row in rows:
        # absolute path
        row['input_pdb'] = os.path.join(params['rawdata_pdb_dir'], row['relative_path'])
        # parse files and sequences
        key = (row['pdb_id'], row['model_no'])
        if key in pdb_chains:
            chain_ids = row['chain_ids'].split(',')
            sub = pdb_chains[key].sort_values('chain_id')
            sub = sub[sub['chain_id'].isin(chain_ids)]
            # list type
            row['chain_seqs'] = sub['chain_seq'].to_list()
            row['chain_pdbs'] = sub['chain_pdb'].str.cat(sep=',')
            row['chain_nos'] = row['chain_combo'].split('-')
        else:
            print(f"###Warning: no match {row}")
        yield row


def build_data(rows, params, meta):
    '''
    build *.csv, *.fasta, *.pdb
    '''
    for row in rows:
        pdb_id = row['pdb_id']
        endpoint = os.path.join(params['table_name'], pdb_id[:2], pdb_id) 
        outdir = os.path.join(params['output_pdb_dir'], endpoint)
        Path(outdir).mkdir(parents=True, exist_ok=True)
        combo_id = row['combo_id']

        # 1. export combo pdb
        pdb_file_name = f"{combo_id}.pdb"
        pdb_file = os.path.join(outdir, pdb_file_name)
        if params['overwrite'] or not os.path.isfile(pdb_file):
            p = Pdb(row['input_pdb'], None, params['verbose'])
            p.load_structure()
            status = p.to_chain_pdb(row['model_no'], row['chain_nos'], pdb_file)
            meta['combo_pdb_succeed'] += 1
        else:
            meta['combo_pdb_skip'] += 1

        # 2. export complex sequence in csv
        csv_file_name = f"{combo_id}.csv"
        csv_file = os.path.join(outdir, csv_file_name)
        if params['overwrite'] or not os.path.isfile(csv_file):
            combo_seq = {combo_id: row['chain_seqs'],}        
            status = ProcessSeq.to_colabfold_csv(combo_seq, csv_file)
            meta['combo_csv_'+ status] += 1
        else:
            meta['combo_csv_skip'] += 1
        
        # 3. export complex sequence in fasta
        fasta_file_name = f"{combo_id}.fasta"
        fasta_file = os.path.join(outdir, fasta_file_name)
        if params['overwrite'] or not os.path.isfile(fasta_file):
            combo_seq = {combo_id: row['chain_seqs'],}        
            status = ProcessSeq.to_colabfold_fasta(combo_seq, fasta_file)
            meta['combo_fasta_'+ status] += 1
        else:
            meta['combo_fasta_skip'] += 1

        # prepare record
        bounded_pdb = os.path.join(endpoint, pdb_file_name) if \
            os.path.isfile(pdb_file) else None
        complex_csv = os.path.join(endpoint, csv_file_name) if \
            os.path.isfile(csv_file) else None
        complex_fasta = os.path.join(endpoint, fasta_file_name) if \
            os.path.isfile(fasta_file) else None
        rec = (pdb_id, combo_id, row['chain_combo'], row['chain_pdbs'], \
            bounded_pdb, complex_csv, complex_fasta)
        if params['verbose']:
            print(params['pdb_group'], rec)
        yield rec

##################################
if __name__ == "__main__":
    params.update({
        'overwrite': False,
        'combo': 3,
        'table_name': 'combo3',
    })
    row_iter = retrieve_data(params, meta)
    record_iter = build_data(row_iter, params, meta)
    bc = BuildComplex(params['verbose'])
    bc.insert_batch_records(record_iter, params['table_name'])

    for col_name in ('combo_pdb', 'combo_csv', 'combo_fasta'):
        validate_data(col_name, params, meta)
    footer(meta)
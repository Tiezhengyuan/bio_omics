'''
requirements:
    - The table meta_pdb is ready.
    - meta_chain_pdb is ready
example:
    - parallel -j8 python src/pipelines/meta_combo2.py ::: {0..117}
    - python src/pipelines/meta_combo2.py 0
function:
    - split PDB of complex to paired chains PDB.
    - export complex csv and fasta for structure predication
note:
    - only consider model.id == 0
'''
from helper import *
from itertools import combinations
from bioomics import ProcessSeq


@header
def retrieve_data(params, meta):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        SELECT C.pdb_id, C.chain_id, C.chain_pdb, I.model_no, I.chain_no,
            P.relative_path, S.chain_seq
        FROM meta_chain_pdb     C
        LEFT JOIN pdb_chainid   I ON C.chain_id=I.chain_id
        LEFT JOIN pdb_chain_seq S ON C.chain_id=S.chain_id
        LEFT JOIN meta_pdb      P ON C.pdb_id=P.pdb_id
        WHERE P.pdb_group = {pdb_group}
            AND P.relative_path IS NOT NULL
            AND S.chain_seq IS NOT NULL
            AND I.model_no = 0
            AND C.pdb_id NOT IN (
                SELECT pdb_id FROM {table_name}
            )
            AND C.pdb_id IN (
                SELECT pdb_id from view_antibody
            )
    ;"""
    df = QueryComplex(params['verbose']).list_data(query, True)
    meta['pdb_chains'] = len(df)

    if len(df) > 0:
        print(f"Filter {len(df)} PDB files for the next data processing...")
        g = df.groupby(['pdb_id', 'relative_path', 'model_no'])
        for (pdb_id, relative_path, model_no), sub in g:
            if len(sub) >= params['combo']:
                pdb_file = os.path.join(params['pdb_dir'], relative_path)
                if os.path.isfile(pdb_file):
                    sub2 = sub.sort_values('chain_no').set_index('chain_no')
                    meta['input_chains'] += 1
                    yield {
                        'input_pdb': pdb_file,
                        'pdb_id': pdb_id,
                        'model_no': model_no,
                        'chain_rows': sub2,
                    }
                else:
                    meta['no_pdb_file'] += 1
            else:
                meta['no_combo'] += 1

def build_data(input_iter, params, meta):
    for input in input_iter:
        pdb_id = input['pdb_id']
        endpoint = os.path.join(params['table_name'], pdb_id[:2], pdb_id) 
        outdir = os.path.join(params['data_dir'], endpoint)
        Dir(outdir).init_dir()

        p = Pdb(input['input_pdb'], None, params['verbose'])
        p.load_structure()

        model_no = input['model_no']
        chain_rows = input['chain_rows']
        ids = list(chain_rows.index)
        for _ids in combinations(ids, params['combo']):
            _ids = list(_ids)
            chain_combo = '-'.join(_ids)
            combo_id = f"{pdb_id}_{model_no}_{chain_combo}"

            # combo pdb
            pdb_file_name = f"{combo_id}.pdb"
            pdb_file = os.path.join(outdir, pdb_file_name)
            if params['overwrite'] or not os.path.isfile(pdb_file):
                status = p.to_chain_pdb(model_no, _ids, pdb_file)
                meta['combo_pdb_succeed'] += 1
            else:
                meta['combo_pdb_skip'] += 1

            # export complex sequences
            combo_seq = {
                combo_id: list(chain_rows.loc[_ids, 'chain_seq']),
            }        
            # complex sequence in csv
            csv_file_name = f"{combo_id}.csv"
            csv_file = os.path.join(outdir, csv_file_name)
            if params['overwrite'] or not os.path.isfile(csv_file):
                status = ProcessSeq.to_colabfold_csv(combo_seq, csv_file)
                meta['combo_csv_'+ status] += 1
            else:
                meta['combo_csv_skip'] += 1
            # complex sequence in fasta
            fasta_file_name = f"{combo_id}.fasta"
            fasta_file = os.path.join(outdir, fasta_file_name)
            if params['overwrite'] or not os.path.isfile(fasta_file):
                status = ProcessSeq.to_colabfold_fasta(combo_seq, fasta_file)
                meta['combo_fasta_'+ status] += 1
            else:
                meta['combo_fasta_skip'] += 1

            # prepare record
            chain_pdbs = chain_rows.loc[_ids, 'chain_pdb'].str.cat(sep=',')
            bounded_pdb = os.path.join(endpoint, pdb_file_name) if \
                os.path.isfile(pdb_file) else None
            complex_csv = os.path.join(endpoint, csv_file_name) if \
                os.path.isfile(csv_file) else None
            complex_fasta = os.path.join(endpoint, fasta_file_name) if \
                os.path.isfile(fasta_file) else None
            rec = (pdb_id, combo_id, chain_combo, chain_pdbs, \
                bounded_pdb, complex_csv, complex_fasta)
            if params['verbose']:
                print(params['pdb_group'], rec)
            yield rec

##################################
if __name__ == "__main__":
    params = {
        'pdb_group': sys.argv[1],
        'verbose': True,
        'overwrite': False,
        'chunk_size': 100,
        'pdb_dir': os.getenv('pdb_dir'),
        'data_dir': os.getenv('data_dir'),
        'table_name': 'meta_combo2',
        'combo': 2,
    }
    input_iter = retrieve_data(params, meta)
    record_iter = build_data(input_iter, params, meta)
    bc = BuildComplex(record_iter, params['verbose'])
    bc.insert_batch_records(params['table_name'])

    for col_name in ('combo_pdb', 'combo_csv', 'combo_fasta'):
        validate_data(col_name, params, meta)
    footer(meta)
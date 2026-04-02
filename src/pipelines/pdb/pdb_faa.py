'''
requirements:
    - PDB data
    - table meta_pdb
example: 
    - parallel -j4 python app.py pdb_faa ::: {0..237}
    - python app.py pdb_faa 0
functions:
    - retrieve pdb to fasta
    - blastp against uniprot_sprot
'''

from src.helper import *
from bioomics import Pdb, ProcessFasta

def retrieve_pdb(params, meta):
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
    meta['rows'] = len(rows)

    print(f"{len(rows)} PDB files ready for processing...")
    for row in rows:
        # load structure
        pdb_id = row['pdb_id']
        pdb_file = os.path.join(params['pdb_dir'], row['relative_path'])
        p = Pdb(pdb_file, None, params['verbose'])
        p.load_structure()

        # retrieve sequences by models
        model_chains = p.parse_model_chains_seq()
        if model_chains:
            # prepare direcotry
            endpoint = os.path.join( params['table_name'], pdb_id[:2], pdb_id)
            outdir = os.path.join(params['data_dir'], endpoint)
            Dir(outdir).init_dir()
            
            for model_no, chains_seq in model_chains.items():
                file_name = f"{pdb_id}_{model_no}.faa"
                outfile = os.path.join(outdir, file_name)
                status = ProcessFasta(outfile).row_seq(chains_seq)
                meta[status] += 1
                if os.path.isfile(outfile):
                    relative_path = os.path.join(endpoint, file_name)
                    yield pdb_id, model_no, relative_path
                else:
                    yield pdb_id, model_no, None
        else:
            meta['no_seq'] += 1
            yield pdb_id, None, None


if __name__ == "__main__":
    params = {
        'pdb_group': sys.argv[1],
        'verbose': True,
        'pdb_dir': os.getenv('pdb_dir'),
        'data_dir': os.getenv('data_dir'),
        'chunk_size': 50,
        'table_name': 'pdb_faa',
    }

    # retrieve pdb_ids
    record_iter = retrieve_pdb(params, meta)
    bc = BuildComplex(record_iter, params['verbose'], params['chunk_size'])
    bc.insert_batch_records(params['table_name'])

    footer(meta)


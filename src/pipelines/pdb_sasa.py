'''
requirments:
    - PDB raw data
    - table meta_combo2_pdb
example: 
    - parallel -j8 python src/pipelines/pdb_sasa.py ::: {0..117}
    - python src/pipelines/pdb_sasa.py 0
functions:
    - calculate SASA and analyze PPI
    - put meta data to table pdb_sasa
'''

from helper import *
from bioomics import SasaProcess

def retrieve_query(pdb_group, args, meta):
    table_name = args['table_name']
    query = f"""
        SELECT S.pdb_id, S.combo_id, S.combo_pdb, S.chain_pdbs
        FROM meta_combo2   S
        LEFT JOIN meta_pdb P ON S.pdb_id = P.pdb_id
        WHERE P.pdb_group = {pdb_group}
            AND P.relative_path IS NOT NULL
            AND S.combo_id IS NOT NULL
            AND S.combo_pdb IS NOT NULL
            AND S.chain_pdbs IS NOT NULL
            AND S.combo_id NOT IN(
                SELECT combo_id FROM {table_name}
            );
    """
    rows = QueryComplex(args['verbose']).list_data(query)
    meta['rows'] = len(rows)
    meta['pdb_group'] = pdb_group
    print(f"Group {pdb_group}: Collect combo2 chains from {len(rows)} PDB ...")
    return rows

def cal_sasa(rows, args, meta):
    for row in rows:
        pdb_id = row['pdb_id']
        endpoint = os.path.join(args['table_name'], pdb_id[:2], pdb_id)
        outdir = os.path.join(args['data_dir'], endpoint)
        Dir(outdir).init_dir()

        file_name =  f"{row['combo_id']}.pkl"
        pfile = os.path.join(outdir, file_name)
        if args['overwrite'] or (not os.path.isfile(pfile)):
            try:
                combo_pdb_file = os.path.join(args['data_dir'], row['combo_pdb'])
                chain_pdbs = row['chain_pdbs'].split(',')
                chain_pdb_files = [os.path.join(args['data_dir'], i) for i in chain_pdbs]
                pdb_files_pool = [(combo_pdb_file, chain_pdb_files),]
                SasaProcess(verbose=args['verbose']).process(pdb_files_pool, pfile)
                meta['sasa_succeed'] += 1
            except Exception as e:
                print(e)
                meta['sasa_fail'] += 1
        else:
            meta['skip'] += 1
        
        # build record
        relative_path = os.path.join(endpoint, file_name) if os.path.isfile(pfile) else None
        rec = (pdb_id, row['combo_id'], relative_path)
        yield rec


if __name__ == "__main__":
    args = {
        'verbose': True,
        'overwrite': False,
        'chunk_size': 20,
        'data_dir': os.getenv('data_dir'),
        'table_name': 'pdb_sasa',
    }
    pdb_group = sys.argv[1]
    rows = retrieve_query(pdb_group, args, meta)
    # calculate delta total SASA
    record_iter = cal_sasa(rows, args, meta)

    print(f"Try to put meta data to table {args['table_name']}")
    bc = BuildComplex(record_iter, args['verbose'], args['chunk_size'])
    bc.insert_batch_records(args['table_name'])

    print(meta)

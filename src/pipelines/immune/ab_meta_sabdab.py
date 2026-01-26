'''
Download data from database SAbDab https://opig.stats.ox.ac.uk
example: 
    - python src/pipelines/ab_sabdab.py
functions:
    - download data from SAbDab
    - put data to table ab_sabdab
'''

from helper import *
from bioomics import Sabdab


if __name__ == '__main__':
    meta = defaultdict(int)
    load_dotenv()
    args = {
        'verbose': True,
        'overwrite': False,
        'chunk_size': 50,
        'download_dir': os.getenv("download_dir"),
        'table_name': 'meta_sabdab',
    }
    sab = Sabdab(args['download_dir'], args['overwrite'])

    # retrieve list of pdb_id from sabdab_pdbs_list.txt
    # Warning: the list pdb in this file doesn't cover all pdb in sabdab.
    sabdab_ids = sab.retrieve_sabdab_pdb_ids()
    meta['sabdab_pdb_ids'] = len(sabdab_ids)

    # filter pdb_ids
    table_name = 'meta_sabdab'
    pdb_ids = QueryComplex(args['verbose']).unique_values(table_name, 'pdb_id')
    if pdb_ids:
        sabdab_ids = list(set(sabdab_ids).difference(pdb_ids))
    meta['filtered_pdb_ids'] = len(sabdab_ids)

    # download summary/*.tsv
    record_iter = sab.download_summary(sabdab_ids, meta)

    # insertion
    bc = BuildComplex(record_iter, args['verbose'], args['chunk_size'])
    bc.empty_table(args['table_name'])
    bc.insert_batch_records(args['table_name'])

    print(meta)

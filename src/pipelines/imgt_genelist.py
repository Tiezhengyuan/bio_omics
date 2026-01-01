'''
requirements:
    - IMGT GeneList
example:
    - python src/pipelines/imgt_genelist.py
functions:
    - build table imgt_genelist
'''

from helper import *
from bioomics import ImgtGenedb

if __name__ == "__main__":
    meta = defaultdict(int)
    load_dotenv()
    args = {
        'verbose':True,
        'imgt_dir': os.getenv("imgt_dir"),
    }

    ig = ImgtGenedb(args['imgt_dir'], args['verbose'])
    genelist = ig.parse_genelist(to_df=True)
    records = genelist.to_dict(orient='records')
    meta['rows'] = len(records)

    table_name = 'imgt_genelist'
    print(f"Try to insert genelist into table {table_name}")
    bc = BuildComplex(records, args['verbose'])
    bc.empty_table(table_name)
    bc.insert_imgt_genelist()

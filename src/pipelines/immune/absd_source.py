'''
example: 
    - parallel -j4 python src/pipelines/absd_source.py ::: {0..4}
    - python src/pipelines/absd_source.py 4
functions:
    - download data from ABSD
    - put data to table absd_source
'''
from dotenv import load_dotenv
from collections import defaultdict
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import BuildComplex, Absd


if __name__ == "__main__":
    meta = defaultdict(int)
    load_dotenv()
    args = {
        'verbose': True,
        'chunk_size': 200,
        'absd_dir': os.getenv('absd_dir'),
    }
    data_dir = os.path.join(args['absd_dir'], 'data')
    a = Absd(data_dir, args['verbose'])
    files = [i[1] for i in a.scan_files('fasta')]
    idx = int(sys.argv[1])
    record_iter = a.scan_data(files[idx])

    # insertion
    table_name = 'absd_source'
    bc = BuildComplex(record_iter, args['verbose'], args['chunk_size'])
    bc.insert_absd_source()
    
    print(meta)


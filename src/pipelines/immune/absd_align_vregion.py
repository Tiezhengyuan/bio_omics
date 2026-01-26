#!/usr/bin/python
'''
requirements:
    - table absd_faa, absd_igblastp
    - igblastp ready
example:
    - parallel -j8 python src/pipelines/absd_align_vregion.py ::: {0..460}
    - python src/pipelines/absd_align_vregion.py 0
function:
    - build table pdb and insert data into table absd_align_vregion
'''

from helper import *

def retrieve_data(idx, args, meta):
    end = (idx + 1 ) * 10_000
    start = end - 10_000
    query = f"""
        SELECT F.seq_id, F.faa_file, I.aln_file 
        FROM absd_faa F, absd_igblastp I
        WHERE F.seq_id >= {start}
            AND F.seq_id < {end}
            AND F.seq_id = I.seq_id
            AND F.faa_file IS NOT NULL
            AND I.aln_file IS NOT NULL
            AND F.seq_id NOT IN(
                SELECT DISTINCT seq_id
                FROM absd_align_vregion
            );
    """
    rows = QueryComplex(args['verbose']).list_data(query)
    meta['idx'] = int(idx)
    meta['rows'] = len(rows)
    return rows

if __name__ == "__main__":
    args = {
        'verbose': True,
    }
    idx = int(sys.argv[1])
    rows = retrieve_data(idx, args, meta)

    table_name = 'absd_align_vregion'
    print(f"Try to insert data into {table_name}")
    bc = BuildComplex(rows, args['verbose'])
    bc.insert_absd_align_vregion()

    print(meta)


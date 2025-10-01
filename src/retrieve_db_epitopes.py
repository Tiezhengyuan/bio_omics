'''
retrieve epitopes from IEDB
'''
import pandas as pd
import json
import gzip
from zipfile import ZipFile
import mysql.connector
import os


def retrieve(outdir:str, condtion:str) -> tuple:
    keys = ['id', 'start', 'end', 'seq']
    n, epitopes = 0, {}

    db_params = {
        'host': "localhost",
        'user': "root",
        'passwd': "root",
        'database': 'protein',
    }
    conn = mysql.connector.connect(**db_params)
    cursor = conn.cursor()
    cursor.execute(f"""
        SELECT a.curated_epitope_id,
            b.starting_position, b.ending_position, b.mol1_seq,
            c.sequence, c.accession
        FROM curated_epitope a, object b, source c
        WHERE a.e_object_id = b.object_id
            AND b.mol2_source_id = c.source_id
            AND b.starting_position is NOT NULL
            AND c.sequence is NOT NULL
            AND LEFT(c.accession, 1) = '{condition}'
    """)
    for row in cursor:
        row = [i if isinstance(i, str) else int(i) for i in row]
        acc = row[-1]
        if acc not in epitopes:
            epitopes[acc] = {
                'pro_seq': row[-2],
                'epitopes': {},
            }
        rec = {a:b for a,b in zip(keys, row[:-2])}
        epitopes[acc]['epitopes'][rec['id']] = rec
        n += 1
        if n % 10_000 == 0:
            print(f"{condition}:{int(n/1000)}", end=', ')
            # break
    conn.close()

    # export
    print('Number of proteins: ', len(epitopes))    
    outfile = os.path.join(outdir, f'epitope_{condition}.json')
    with open(outfile, 'w') as f:
        json.dump(epitopes, f, indent=4, sort_keys=True)
    return outfile, n, len(epitopes)

num_epitopes, num_proteins = 0, 0
outfiles = []
outdir = '/home/yuan/data/omics_data/epitope/mysql'
conditions = [chr(i) for i in range(65, 91)] + list(range(0, 10))
for condition in conditions:
    outfile, n, p = retrieve(outdir, condition)
    outfiles.append(outfile)
    num_epitopes += n
    num_proteins += p
print('Total number of epitopes: ', num_epitopes)
print('Total number of proteins: ', num_proteins)


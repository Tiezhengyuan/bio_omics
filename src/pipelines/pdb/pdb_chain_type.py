#!/usr/bin/python
'''
requirements:
    - table meta_pdb
example:
    - parallel -j8 bomics pdb_chain_type ::: {0..241}
    - bomics pdb_chain_type 0
function:
    - identify chain type according pdb file
'''

from src.helper import *
from bioomics import Pdb

@header
def retrieve_pdb(params, meta):
    qc = QueryComplex(params['verbose'])
    rows = qc.meta_pdb(params['pdb_group'], params['table_name'])
    print(f"{len(rows)} PDB files ready...")
    meta['rows'] += len(rows)
    for row in rows:
        pdb_file = params['rawdata_dir'] / row['relative_pdb']
        if pdb_file.is_file():
            row['pdb_file'] = pdb_file
            yield row

def detect_chain_type(rows, params, meta):
    for row in rows:
        pdb_id = row['pdb_id']
        meta['pdb_record'] += 1
        p = Pdb(row['pdb_file'], None, params['verbose'])
        p.load_structure()

        chains = p.parse_chains()
        if chains:
            for item in chains:
                counts = p.residue_type(item['chain'])
                rec = (
                    pdb_id,
                    item['chain_id'],
                    counts['type'],
                    counts['total_res'],
                    counts['non_water_res'],
                    counts['aa_res'],
                    counts['dna_res'],
                    counts['rna_res'],
                    counts['water_res'],
                    counts['ion_res'],
                    counts['hetero_res'],
                    counts['other_res']
                )
                meta['chain'] += 1
                yield rec
        else:
            yield (pdb_id,) + (None,) * 11
            meta['null'] += 1


if __name__ == "__main__":
    params.update({
        'chunk_size': 50,
        'table_name': 'pdb_chain_type',
    })
    # retrieve chain types
    rows = retrieve_pdb(params, meta)
    record_iter = detect_chain_type(rows, params, meta)

    bc = BuildComplex(params['verbose'], params['chunk_size'])
    bc.insert_batch_records(record_iter, params['table_name'])

    footer(meta)
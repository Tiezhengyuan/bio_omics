'''
requirements:
    - PDB: meta_pdb, pdb_chain
    - blastp of pdb chain and table align_uniprot
    - alphafold prediction in PDB
example: 
    - parallel -j4 python src/pipelines/uniprot_alphafold_rmsd.py ::: {0..46}
    - python src/pipelines/uniprot_alphafold_rmsd.py 46
functions:
    - similarity between model and sample structure
'''
from collections import defaultdict
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import QueryComplex, AnalyzeModel, BuildComplex


def search_pdb(pdb_group, args, meta):
    query = f"""
        SELECT C.chain_id, C.pdb_id, U.uniprot_acc, 
            C.pdb_file as native_pdb, A.pdb_path as model_pdb
        FROM view_pdb_chain C
        LEFT JOIN meta_pdb M
            ON C.pdb_id = M.pdb_id
        LEFT JOIN align_uniprot U
            ON C.chain_id = U.chain_id
        LEFT JOIN uniprot_meta_alphafold A
            ON U.uniprot_acc = A.uniprot_acc
        WHERE M.pdb_group = {pdb_group}
            AND U.bit_score >= 50
            AND U.evalue <= 1e-10
            AND U.identity >= .3
            AND A.pdb_path IS NOT NULL
            AND C.pdb_file IS NOT NULL
            AND C.chain_id NOT IN (
                SELECT DISTINCT chain_id
                FROM uniprot_alphafold_rmsd
            );
    """
    rows = QueryComplex(args['verbose']).list_data(query)
    meta['pdb_group'] = int(pdb_group)
    meta['rows'] = len(rows)
    return rows

def calculate_rmsd(rows, meta):
    for row in rows:
        rmsd = None
        try:
            rmsd = AnalyzeModel.rmsd(row['native_pdb'], row['model_pdb'])
            if rmsd:
                meta['succeed_rmsd'] += 1
            else:
                meta['null_rmsd'] += 1    
        except Exception as e:
            meta['fail_rmsd'] += 1
        # wrap record
        rec = (row['chain_id'], row['uniprot_acc'], rmsd,\
            row['native_pdb'], row['model_pdb'])
        yield rec


if __name__ == "__main__":
    args = {
        'verbose': True,
        'chunk_size': 10,
    }
    meta = defaultdict(int)

    # retrieve input
    pdb_group = sys.argv[1]
    rows = search_pdb(pdb_group, args, meta)

    # RMSD
    record_iter = calculate_rmsd(rows, meta)

    # insert data
    bc = BuildComplex(record_iter, args['verbose'], args['chunk_size'])
    bc.insert_uniprot_alphafold_rmsd()

    print(meta)

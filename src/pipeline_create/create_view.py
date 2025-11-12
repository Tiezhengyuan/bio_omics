#!/usr/bin/python
'''
requirements:
    - a empty database is created.
example: python src/pipelines/complex_create.py
functions:
    - build database known as complex2
    - build tables
'''
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)


from bioomics import BuildComplex

SQL_TEXT = """
DROP VIEW IF EXISTS view_pdb;
CREATE VIEW view_pdb AS
    SELECT P.pdb_id, P.pdb_url, B.pdb_file, P.release_date,
        P.name, M.structure_method, M.resolution
    FROM pdb P, meta_pdb B, pdb_method M
    WHERE P.pdb_id = B.pdb_id
        AND P.pdb_id = M.pdb_id
        AND M.resolution IS NOT NULL
;
//
DROP VIEW IF EXISTS view_pdb_chain;
CREATE VIEW view_pdb_chain AS
    SELECT C.chain_id, C.pdb_id, C.model_no, C.chain_no,
        M.structure_method, M.resolution,
        S.chain_seq, F.faa_file, P.pdb_file, D.pdb_group
    FROM pdb_chainid C
    LEFT JOIN pdb_chain_seq S
        ON C.chain_id = S.chain_id
    LEFT JOIN meta_chain_faa F
        ON C.chain_id = F.chain_id
    LEFT JOIN meta_chain_pdb P
        ON C.chain_id = P.chain_id
    LEFT JOIN pdb_method M
        ON C.pdb_id = M.pdb_id
    LEFT JOIN meta_pdb D
        ON C.pdb_id = D.pdb_id
;
//
DROP VIEW IF EXISTS view_antibody;
CREATE VIEW view_antibody AS
    SELECT I.allele_name, I.gene_name, I.gene_family,
        I.chain_type, C.pdb_id, C.chain_id, C.model_no,
        C.chain_no, C.chain_seq, C.structure_method,
        C.resolution, C.faa_file, C.pdb_file, C.pdb_group
    FROM view_pdb_chain C, align_imgt I
    WHERE C.chain_id = I.chain_id
        AND I.evalue <= 1e-10
        AND I.bit_score >= 50
        AND I.isotype = 'IG'
;
//
DROP VIEW IF EXISTS view_pdb_uniprot;
CREATE VIEW view_pdb_uniprot AS
    SELECT A.chain_id, A.pdb_id, A.uniprot_acc, U.name,
        C.chain_no, S.chain_seq, A.identity, U.pro_seq
    FROM align_uniprot A
    LEFT JOIN uniprot U
        ON A.uniprot_acc = U.uniprot_acc
    LEFT JOIN pdb_chainid C
        ON A.chain_id = C.chain_id
    LEFT JOIN pdb_chain_seq S
        ON C.chain_id = S.chain_id
    WHERE A.bit_score >= 50
        AND A.evalue <= 1e-10
        AND A.identity >= .3
        AND C.model_no = 0
        AND S.chain_seq IS NOT NULL
;
//
DROP VIEW IF EXISTS view_antibody_complex;
CREATE VIEW view_antibody_complex AS
    SELECT * FROM view_pdb_uniprot
    WHERE pdb_id IN(
        SELECT DISTINCT pdb_id
        FROM view_antibody
    )
    ORDER BY pdb_id
;
"""


#####################################################################################

if __name__ == "__main__":
    bc = BuildComplex(verbose=True)
    for query in SQL_TEXT.split('//'): 
        res = bc.create_table(query)
        if res is None:
            print(f"ERROR. Check the above SQL: {query}\n\n")

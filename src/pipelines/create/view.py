#!/usr/bin/python
'''
requirements:
    - a empty database is created.
example:
    - python src/pipeline_create/view.py
functions:
    - create views for database COMPLEX2
'''
from create_helper import *

'''
//
DROP VIEW IF EXISTS view_chain_seq;
CREATE VIEW view_chain_seq AS
    WITH first_chain_seq AS (
        SELECT MIN(chain_id) AS first_chain_id, COUNT(*) AS num_seq,
            chain_seq AS unique_seq
        FROM pdb_chain_seq
        GROUP BY chain_seq
        HAVING chain_seq IS NOT NULL
    )
    SELECT C.pdb_id, C.chain_id, F.first_chain_id,
        F.num_seq, C.pro_len, C.chain_seq
    FROM pdb_chain_seq C
    LEFT JOIN first_chain_seq F
    ON C.chain_seq = F.unique_seq
;
//
DROP VIEW IF EXISTS view_protein;
CREATE VIEW view_protein AS
    SELECT S.source_name, S.source_id, P.pro_id, P.pro_seq
    FROM protein_source   S
    LEFT JOIN protein_seq P ON S.pro_id = P.pro_id
;

'''


SQL_TEXT = """
DROP VIEW IF EXISTS view_pdb;
CREATE VIEW view_pdb AS
    SELECT P.pdb_id, P.pdb_url, B.relative_path, B.pdb_group,
        P.release_date, P.name, M.structure_method, M.resolution
    FROM pdb P, meta_pdb B, pdb_method M
    WHERE P.pdb_id = B.pdb_id
        AND P.pdb_id = M.pdb_id
        AND B.relative_path IS NOT NULL
        AND M.resolution IS NOT NULL
;
//
DROP VIEW IF EXISTS view_pdb_chain;
CREATE VIEW view_pdb_chain AS
    SELECT C.chain_id, C.pdb_id, C.model_no, C.chain_no,
        S.first_chain_id, X.proseq_name,X.seq_id,
        D.structure_method, D.resolution, B.avg_bfactor,
        D.release_date, D.pdb_group,
        S.chain_seq, S.pro_len,
        F.relative_path AS relative_faa, P.chain_pdb
    FROM pdb_chainid        C
    LEFT JOIN pdb_chain_seq S ON C.chain_id = S.chain_id
    LEFT JOIN chain_faa     F ON C.chain_id = F.chain_id
    LEFT JOIN chain_pdb     P ON C.chain_id = P.chain_id
    LEFT JOIN view_pdb      D ON C.pdb_id = D.pdb_id
    LEFT JOIN pdb_pro_seq   X ON C.chain_id = X.chain_id
    LEFT JOIN pdb_bfactor   B ON C.pdb_id = B.pdb_id 
        AND C.model_no = B.model_no AND C.chain_no = B.chain_no
    WHERE S.chain_seq IS NOT NULL
;
//
DROP VIEW IF EXISTS view_pdb_uniprot;
CREATE VIEW view_pdb_uniprot AS
    SELECT A.chain_id, A.pdb_id, A.uniprot_acc, U.name,
        C.chain_no, C.chain_seq,
        A.identity, U.pro_seq,
        C.structure_method, C.resolution,
        C.relative_faa, C.chain_pdb
    FROM align_uniprot       A
    LEFT JOIN uniprot        U ON A.uniprot_acc = U.uniprot_acc
    LEFT JOIN view_pdb_chain C ON A.chain_id = C.chain_id
    WHERE A.bit_score >= 50
        AND A.evalue <= 1e-10
        AND A.identity >= .3
        AND C.model_no = 0
;
//
/* human transmembrane proteins */
DROP VIEW IF EXISTS view_htrm;
CREATE VIEW view_htrm AS
    SELECT C.chain_id, C.pdb_id,  C.uniprot_acc,
        C.chain_no, C.name, C.chain_seq, C.pro_seq
    FROM view_pdb_uniprot C
    LEFT JOIN uniprot U
    ON C.uniprot_acc = U.uniprot_acc
    WHERE U.organism LIKE '%Homo sapiens%'
        AND C.uniprot_acc IN (
            SELECT distinct uniprot_acc
            FROM uniprot_keyword
            WHERE keyword LIKE 'Transmembrane'
        )
;
//
/* drug protein */
DROP VIEW IF EXISTS view_drug;
CREATE VIEW view_drug AS
    SELECT C.chain_id, C.pdb_id,  C.uniprot_acc,
        C.chain_no, C.name, D.db_acc AS drug_acc,
        C.chain_seq, C.pro_seq, C.relative_faa, C.chain_pdb
    FROM view_pdb_uniprot     C
    LEFT JOIN uniprot U ON    C.uniprot_acc = U.uniprot_acc
    LEFT JOIN uniprot_dbxrefs D ON C.uniprot_acc = D.uniprot_acc
    WHERE U.organism LIKE '%Homo sapiens%'
        AND D.db_name = 'DrugBank'
;
//
DROP VIEW IF EXISTS view_complex;
CREATE VIEW view_complex AS
    SELECT pdb_id, count(*) AS num_chains, GROUP_CONCAT(DISTINCT 
        chain_type ORDER BY chain_type) AS complex_type
    FROM pdb_chain_type
    GROUP BY pdb_id
    ORDER BY num_chains, complex_type
;
//
DROP VIEW IF EXISTS view_proseq;
CREATE VIEW view_proseq AS
    SELECT 's400' AS proseq_name, seq_id, source, seq_head, seq_len, seq, relative_faa
    FROM proseq_400
    UNION
    SELECT 's1000' AS proseq_name, seq_id, source, seq_head, seq_len, seq, relative_faa
    FROM proseq_1000
    UNION
    SELECT 'long' AS proseq_name, seq_id, source, seq_head, seq_len, seq, relative_faa
    FROM proseq_long;
"""


#####################################################################################

if __name__ == "__main__":
    print('Load env from ', env_path)
    bc = BuildComplex(verbose=True)
    for query in SQL_TEXT.split('//'): 
        res = bc.create_table(query)
        if res is None:
            print(f"ERROR. Check the above SQL: {query}\n\n")

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
    SELECT P.pdb_id, P.pdb_url, B.pdb_file, B.pdb_group,
        P.release_date, P.name, M.structure_method, M.resolution
    FROM pdb P, meta_pdb B, pdb_method M
    WHERE P.pdb_id = B.pdb_id
        AND P.pdb_id = M.pdb_id
        AND M.resolution IS NOT NULL
;
//
DROP VIEW IF EXISTS view_protein;
CREATE VIEW view_protein AS
    SELECT S.source_name, S.source_id, P.pro_id, P.pro_seq
    FROM protein_source S
    LEFT JOIN protein_seq P
    ON S.pro_id = P.pro_id
;
//
DROP VIEW IF EXISTS view_pdb_chain;
CREATE VIEW view_pdb_chain AS
    SELECT C.chain_id, X.pro_id, C.pdb_id, C.model_no,
        C.chain_no, D.structure_method, D.resolution,
        S.chain_seq, F.faa_file, P.pdb_file, D.pdb_group
    FROM pdb_chainid C
    LEFT JOIN pdb_chain_seq S
        ON C.chain_id = S.chain_id
    LEFT JOIN meta_chain_faa F
        ON C.chain_id = F.chain_id
    LEFT JOIN meta_chain_pdb P
        ON C.chain_id = P.chain_id
    LEFT JOIN view_pdb D
        ON C.pdb_id = D.pdb_id
    LEFT JOIN view_protein X
        ON C.chain_id = X.source_id
    WHERE S.chain_seq IS NOT NULL
        AND X.source_name = 'pdb_chain_seq'
;
//
DROP VIEW IF EXISTS view_pdb_uniprot;
CREATE VIEW view_pdb_uniprot AS
    SELECT A.chain_id, A.pdb_id, A.uniprot_acc, U.name,
        C.chain_no, C.chain_seq,
        A.identity, U.pro_seq,
        C.structure_method, C.resolution,
        C.faa_file, C.pdb_file
    FROM align_uniprot A
    LEFT JOIN uniprot U
        ON A.uniprot_acc = U.uniprot_acc
    LEFT JOIN view_pdb_chain C
        ON A.chain_id = C.chain_id
    WHERE A.bit_score >= 50
        AND A.evalue <= 1e-10
        AND A.identity >= .3
        AND C.model_no = 0
;
//
DROP VIEW IF EXISTS view_antibody;
CREATE VIEW view_antibody AS
    SELECT I.allele_name, I.gene_name, I.gene_family,
        I.chain_type, C.pdb_id, C.chain_id, C.model_no,
        C.chain_no, C.pro_id, C.chain_seq,
        C.structure_method, C.resolution,
        C.faa_file, C.pdb_file, C.pdb_group
    FROM view_pdb_chain C, align_imgt I
    WHERE C.chain_id = I.chain_id
        AND I.evalue <= 1e-10
        AND I.bit_score >= 50
        AND I.isotype = 'IG'
;
//
DROP VIEW IF EXISTS view_tcr;
CREATE VIEW view_tcr AS
    SELECT I.allele_name, I.gene_name, I.gene_family,
        I.chain_type, C.pdb_id, C.chain_id, C.model_no,
        C.chain_no, C.pro_id, C.chain_seq, C.structure_method,
        C.resolution, C.faa_file, C.pdb_file, C.pdb_group
    FROM view_pdb_chain C, align_imgt I
    WHERE C.chain_id = I.chain_id
        AND I.evalue <= 1e-10
        AND I.bit_score >= 50
        AND I.isotype = 'TR'
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
//
/* antibody sequences from ABSD aligned with IMGT V-region */
DROP VIEW IF EXISTS view_abseq_imgt;
CREATE VIEW view_abseq_imgt AS
    SELECT S.seq_id, S.specie, I.allele_name, I.gene_name, I.gene_family,
        I.chain_type, I.isotype, S.pro_seq, F.faa_file, A.aln_file
    FROM absd_seq S, absd_faa F, absd_igblastp A, absd_align_imgt I
    WHERE S.seq_id = F.seq_id
        AND S.seq_id = A.seq_id
        AND S.seq_id = I.seq_id
        AND I.evalue <= 1e-10
        AND I.bit_score >= 50
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
    );
//
/* drug protein */
DROP VIEW IF EXISTS view_drug;
CREATE VIEW view_drug AS
    SELECT C.chain_id, C.pdb_id,  C.uniprot_acc,
        C.chain_no, C.name, D.db_acc AS drug_acc,
        C.chain_seq, C.pro_seq, C.faa_file, C.pdb_file
    FROM view_pdb_uniprot C
    LEFT JOIN uniprot U
    ON C.uniprot_acc = U.uniprot_acc
    LEFT JOIN uniprot_dbxrefs D
    ON C.uniprot_acc = D.uniprot_acc
    WHERE U.organism LIKE '%Homo sapiens%'
        AND D.db_name = 'DrugBank'
;
//
/* allergen */
DROP VIEW IF EXISTS view_allergen;
CREATE VIEW view_allergen AS
    SELECT C.chain_id, C.pdb_id,  C.uniprot_acc,
        C.chain_no, C.name, D.db_acc AS allergen_acc,
        C.chain_seq, C.pro_seq, C.faa_file, C.pdb_file
    FROM view_pdb_uniprot C
    LEFT JOIN uniprot U
    ON C.uniprot_acc = U.uniprot_acc
    LEFT JOIN uniprot_dbxrefs D
    ON C.uniprot_acc = D.uniprot_acc
    WHERE U.organism LIKE '%Homo sapiens%'
        AND D.db_name = 'allergome'
;
"""


#####################################################################################

if __name__ == "__main__":
    bc = BuildComplex(verbose=True)
    for query in SQL_TEXT.split('//'): 
        res = bc.create_table(query)
        if res is None:
            print(f"ERROR. Check the above SQL: {query}\n\n")

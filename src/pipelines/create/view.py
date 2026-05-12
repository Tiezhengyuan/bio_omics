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


SQL_TEXT = """
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
/* allergen */
DROP VIEW IF EXISTS view_allergen;
CREATE VIEW view_allergen AS
    SELECT C.chain_id, C.pdb_id,  C.uniprot_acc,
        C.chain_no, C.name, D.db_acc AS allergen_acc,
        C.chain_seq, C.pro_seq, C.relative_faa , C.relative_pdb
    FROM view_pdb_uniprot     C
    LEFT JOIN uniprot         U ON C.uniprot_acc = U.uniprot_acc
    LEFT JOIN uniprot_dbxrefs D ON C.uniprot_acc = D.uniprot_acc
    WHERE D.db_name = 'allergome'
        AND U.organism LIKE '%Homo sapiens%'
;
//
DROP VIEW IF EXISTS view_tcr;
CREATE VIEW view_tcr AS
    SELECT I.allele_name, I.gene_name, I.gene_family,
        I.chain_type, C.pdb_id, C.chain_id, C.model_no,
        C.chain_no, C.pro_id, C.chain_seq, C.structure_method,
        C.resolution, C.relative_faa, C.relative_pdb, C.pdb_group
    FROM view_pdb_chain C, align_vregion I
    WHERE C.chain_id = I.chain_id
        AND I.evalue <= 1e-10
        AND I.bit_score >= 50
        AND I.isotype = 'TR'
;       
//
DROP VIEW IF EXISTS view_ig;
CREATE VIEW view_ig AS
    SELECT I.allele_name, I.specie, I.gene_name, I.gene_family,
        I.chain_type, C.pdb_id, C.chain_id, C.first_chain_id,
        C.model_no, C.chain_no, C.pro_id, C.chain_seq, C.pro_len,
        C.release_date, C.structure_method, C.resolution, C.avg_bfactor,
        C.relative_faa, C.chain_pdb, C.pdb_group
    FROM view_pdb_chain     C
    LEFT JOIN align_vregion I ON C.chain_id = I.chain_id
    WHERE C.model_no = 0
        AND I.evalue <= 1e-5
        AND I.bit_score >= 50
        AND I.isotype = 'IG'
;
//
//
DROP VIEW IF EXISTS view_ab_dssp;
CREATE VIEW view_ab_dssp AS
    WITH tmp AS (
        SELECT I.allele_name, I.specie, I.gene_name, I.gene_family,
            I.chain_type, C.pdb_id, C.chain_id, C.model_no, C.chain_no,
            S.chain_seq, S.pro_len
        FROM align_vregion I, pdb_chainid C, pdb_chain_seq S
        WHERE I.chain_id = C.chain_id
            AND I.chain_id = S.chain_id
            AND I.evalue <= 1e-10
            AND I.bit_score >= 100
            AND I.isotype = 'IG'
    )
    SELECT A.pdb_id, A.chain_id, A.allele_name, A.specie,
        A.gene_name, A.gene_family, A.chain_type, B.avg_bfactor,
        A.chain_seq, B.bfactor, D.dssp_seq
    FROM tmp A
    LEFT JOIN pdb_bfactor B ON A.pdb_id=B.pdb_id 
        AND A.model_no=B.model_no AND A.chain_no=B.chain_no
    LEFT JOIN pdb_dssp D ON A.pdb_id=D.pdb_id
        AND A.model_no=D.model_no AND A.chain_no=D.chain_no
    WHERE B.bfactor is not null
        AND D.dssp_seq is not null
;
//
DROP VIEW IF EXISTS view_antibody_complex;
CREATE VIEW view_antibody_complex AS
    SELECT * FROM view_pdb_uniprot
    WHERE pdb_id IN (
        SELECT DISTINCT pdb_id FROM view_antibody
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
"""


#####################################################################################

if __name__ == "__main__":
    print('Load env from ', env_path)
    bc = BuildComplex(verbose=True)
    for query in SQL_TEXT.split('//'): 
        res = bc.create_table(query)
        if res is None:
            print(f"ERROR. Check the above SQL: {query}\n\n")

#!/usr/bin/python
'''
requirements:
    - a empty database is created.
example:
    - python src/pipeline_create/view_immune.py
functions:
    - create views for database COMPLEX2
'''
from create_helper import *

'''
//
/* allergen */
DROP VIEW IF EXISTS view_allergen;
CREATE VIEW view_allergen AS
    SELECT C.chain_id, C.pdb_id,  C.uniprot_acc,
        C.chain_no, C.name, D.db_acc AS allergen_acc,
        C.chain_seq, C.pro_seq, C.relative_faa , C.relative_pdb
    FROM view_pdb_uniprot C
    LEFT JOIN uniprot U
    ON C.uniprot_acc = U.uniprot_acc
    LEFT JOIN uniprot_dbxrefs D
    ON C.uniprot_acc = D.uniprot_acc
    WHERE D.db_name = 'allergome'
        AND U.organism LIKE '%Homo sapiens%'
        
;
'''


SQL_TEXT = """
DROP VIEW IF EXISTS view_antibody;
CREATE VIEW view_antibody AS
    SELECT I.allele_name, I.gene_name, I.gene_family,
        I.chain_type, C.pdb_id, C.chain_id, C.model_no,
        C.chain_no, C.pro_id, C.chain_seq, C.release_date,
        C.structure_method, C.resolution,
        C.relative_faa, C.relative_pdb, C.pdb_group
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
        C.resolution, C.relative_faa, C.relative_pdb, C.pdb_group
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
    WHERE pdb_id IN (
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
DROP VIEW IF EXISTS view_inn_pdb;
CREATE VIEW view_inn_pdb AS 
    SELECT B.inn_chain_id, B.pdb_chain_id, P.pdb_id,
        I.inn_number, A.common_name, A.inn_name,
        A.receptor_type, A.receptor_description, I.chain_seq
    FROM imgt_inn_blastp B
    LEFT JOIN pdb_chainid P
        ON B.pdb_chain_id = P.chain_id
    LEFT JOIN imgt_inn_chain I
        ON B.inn_chain_id = I.chain_id
    LEFT JOIN imgt_inn A
        ON I.inn_number = A.inn_number
    WHERE B.identity = 1
;
"""


#####################################################################################

if __name__ == "__main__":
    print('Load env from ', env_path)
    for query in SQL_TEXT.split('//')[:1]:
        bc = BuildComplex(verbose=True)
        res = bc.execute_sql(query)
        if res is None:
            print(f"ERROR. Check the above SQL: {query}\n\n")
            break
#!/usr/bin/python
'''
requirements:
    - a empty database is created.
example: python src/pipelines/create_imgt.py
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
CREATE TABLE imgt_genelist (
    species             VARCHAR(50),
    gene_db             VARCHAR(20),
    gene_function       VARCHAR(20),
    gene_definition     VARCHAR(200),
    number_alleles      VARCHAR(5),
    chromosome          VARCHAR(5),
    chromosomal_loc     VARCHAR(20),
    ligm_db_ref_gene    VARCHAR(100),
    ncbi_gene           VARCHAR(20),
    chain               VARCHAR(10),
    domain              VARCHAR(5)
);
//
DROP TABLE IF EXISTS imgt_gene;
CREATE TABLE imgt_gene (
    record_id       VARCHAR(100),
    specie          VARCHAR(100),
    isotype         VARCHAR(10),
    chain_type      VARCHAR(10),
    gene_name       VARCHAR(20),
    allele_name     VARCHAR(20),
    region_name     VARCHAR(50)
);
//
DROP TABLE IF EXISTS imgt_geneseq;
CREATE TABLE imgt_geneseq (
    record_id       VARCHAR(100),
    seq_type        VARCHAR(100),
    seq             TEXT
);
//
DROP TABLE IF EXISTS msa_vregion;
CREATE TABLE msa_vregion (
    chain_id        VARCHAR(20),
    chain_types     VARCHAR(20),
    query_fa        VARCHAR(100),
    aln_file        VARCHAR(100),
    pickle_file     VARCHAR(100),
    UNIQUE (chain_id, chain_types)
);
//
DROP TABLE IF EXISTS imgt_inn;
CREATE TABLE imgt_inn (
    inn_number              VARCHAR(50),
    cas_number              VARCHAR(100),
    common_name             VARCHAR(100),
    inn_name                VARCHAR(100),
    proposed_list           VARCHAR(50),
    recommended_list        VARCHAR(50),
    receptor_type	        VARCHAR(50),
    receptor_description    VARCHAR(100),
    species    	            VARCHAR(50),
    molecular_formula       VARCHAR(50),
    glycosylation_sites     VARCHAR(200),
    inn_file                VARCHAR(100),
    CONSTRAINT pk_inn PRIMARY KEY (inn_number)
);
//
DROP TABLE IF EXISTS imgt_inn_chain;
CREATE TABLE imgt_inn_chain (
    inn_number              VARCHAR(50),
    chain_id                VARCHAR(50),
    chain_description       VARCHAR(200),
    pro_len                 INT,
    chain_seq               TEXT,
    CONSTRAINT pk_inn_chain PRIMARY KEY (chain_id),
    CONSTRAINT fk_chain_inn FOREIGN KEY (inn_number)
    REFERENCES imgt_inn(inn_number) ON DELETE CASCADE
);
//
DROP TABLE IF EXISTS imgt_inn_region;
CREATE TABLE imgt_inn_region (
    chain_id        VARCHAR(50),
    region_name     VARCHAR(50),
    region_type     VARCHAR(50),
    start           INT,
    end             INT,
    seq             TEXT
);
//
DROP TABLE IF EXISTS imgt_inn_cdomain;
CREATE TABLE imgt_inn_cdomain (
    chain_id        VARCHAR(50),
    domain_name     VARCHAR(50),
    gene_allele     TEXT,
    seq             TEXT,
    seq_align       TEXT
);
//
DROP TABLE IF EXISTS imgt_inn_vdomain;
CREATE TABLE imgt_inn_vdomain (
    chain_id            VARCHAR(50),
    domain_name         VARCHAR(50),
    cdr_length          VARCHAR(50),
    fr_length           VARCHAR(50),
    gene_allele         TEXT,
    seq                 TEXT,
    seq_cdr             TEXT,
    seq_align           TEXT,
    seq_align_cdr       TEXT,
    seq_cdr_mask        TEXT,
    seq_align_cdr_mask  TEXT
);
//
DROP TABLE IF EXISTS imgt_inn_cdr;
CREATE TABLE imgt_inn_cdr (
    chain_id        VARCHAR(50),
    cdr_type        VARCHAR(50),
    seq             VARCHAR(200),
    seq_align       VARCHAR(200),
    align_start     INT,
    align_end       INT,
    align_length    INT
);
"""

#####################################################################################

if __name__ == "__main__":
    bc = BuildComplex(verbose=True)
    for query in SQL_TEXT.split('//'): 
        res = bc.create_table(query)
        if res is None:
            print(f"ERROR. Check the above SQL: {query}\n\n")

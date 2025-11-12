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
    germline_id         INT PRIMARY KEY,
    specie              VARCHAR(100),
    isotype             VARCHAR(10),
    chain_type          VARCHAR(10),
    gene_name           VARCHAR(20),
    allele_name         VARCHAR(20),
    region_name         VARCHAR(50)
);
//
DROP TABLE IF EXISTS imgt_geneseq;
CREATE TABLE imgt_geneseq (
    germline_id     INT,
    seq_type        VARCHAR(100),
    seq             TEXT,
    UNIQUE (germline_id, seq_type),
    CONSTRAINT fk_imgt_genseq_gene FOREIGN KEY (germline_id)
    REFERENCES imgt_gene(germline_id) ON DELETE CASCADE
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
"""

#####################################################################################

if __name__ == "__main__":
    bc = BuildComplex(verbose=True)
    for query in SQL_TEXT.split('//'): 
        res = bc.create_table(query)
        if res is None:
            print(f"ERROR. Check the above SQL: {query}\n\n")

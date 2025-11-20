#!/usr/bin/python
'''
requirements:
    - an empty database is created.
example: python src/pipelines/create_absd.py
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
DROP TABLE IF EXISTS absd_source;
CREATE TABLE absd_source (
    absd_id         VARCHAR(100),
    dxref_id        VARCHAR(200),
    source          VARCHAR(20),
    vgene_family    VARCHAR(20),
    record_id       TEXT,
    header          TEXT
);
DROP TABLE IF EXISTS absd_seq;
CREATE TABLE absd_seq (
    seq_id          INT PRIMARY KEY AUTO_INCREMENT,
    absd_id         VARCHAR(100),
    specie          VARCHAR(100),
    version         VARCHAR(20),
    pro_len         INT,
    pro_seq         TEXT NOT NULL,
    record_id       TEXT
);
//
DROP TABLE IF EXISTS absd_faa;
CREATE TABLE absd_faa (
    seq_id         INT,
    faa_file       VARCHAR(100),
    CONSTRAINT pk_absd_faa PRIMARY KEY (seq_id),
    CONSTRAINT fk_faa_seq FOREIGN KEY (seq_id)
        REFERENCES ab_absd_seq(seq_id)
);
DROP TABLE IF EXISTS absd_igblastp;
CREATE TABLE absd_igblastp (
    seq_id         INT,
    aln_file       VARCHAR(100),
    CONSTRAINT pk_absd_igblastp PRIMARY KEY (seq_id),
    CONSTRAINT fk_igblastp_seq FOREIGN KEY (seq_id)
        REFERENCES ab_absd_seq(seq_id)
);
//
DROP TABLE IF EXISTS absd_align_imgt;
CREATE TABLE absd_align_imgt(
    seq_id        VARCHAR(20),
    allele_name     VARCHAR(20),
    gene_name       VARCHAR(20),
    gene_family     VARCHAR(20),
    chain_type      VARCHAR(10),
    isotype         VARCHAR(10),
    evalue          DOUBLE,
    bit_score       INT,
    CONSTRAINT pk_absd_imgt PRIMARY KEY (seq_id)
);
//
DROP TABLE IF EXISTS absd_align_vregion;
CREATE TABLE absd_align_vregion (
    seq_id        VARCHAR(20),
    region_name     VARCHAR(20),
    identity        FLOAT,
    start           INT,
    end             INT,
    seq             TEXT,
    CONSTRAINT pk_absd_vregion PRIMARY KEY (seq_id, region_name)
);
"""

#####################################################################################

if __name__ == "__main__":
    bc = BuildComplex(verbose=True)
    for query in SQL_TEXT.split('//'): 
        res = bc.create_table(query)
        if res is None:
            print(f"ERROR. Check the above SQL: {query}\n\n")

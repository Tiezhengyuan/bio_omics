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
CREATE TABLE pdb(
    pdb_id          VARCHAR(10),
    pdb_url         VARCHAR(50) NOT NULL,
    release_date    DATE NOT NULL,
    name            VARCHAR(500) NOT NULL,
    CONSTRAINT pk_pdb PRIMARY KEY (pdb_id)
);
//
DROP TABLE IF EXISTS meta_pdb
CREATE TABLE meta_pdb(
    pdb_id          VARCHAR(10),
    pdb_file        VARCHAR(60) NOT NULL,
    pdb_group       INT,
    UNIQUE (pdb_id),
    CONSTRAINT fk_meta_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id)  ON DELETE CASCADE
);
//
CREATE TABLE pdb_compound(
    pdb_id          VARCHAR(10),
    compound_no     INT,
    chains          VARCHAR(200),
    compound_text   TEXT,
    CONSTRAINT pk_compound PRIMARY KEY (pdb_id, compound_no),
    CONSTRAINT fk_compound_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id)  ON DELETE CASCADE
);
//
CREATE TABLE pdb_method(
    pdb_id              VARCHAR(10),
    structure_method    VARCHAR(100),
    resolution          FLOAT,
    CONSTRAINT pk_method PRIMARY KEY (pdb_id),
    CONSTRAINT fk_method_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id)  ON DELETE CASCADE
);
//
CREATE TABLE pdb_keyword(
    pdb_id          VARCHAR(10),
    keyword         VARCHAR(200),
    UNIQUE (pdb_id, keyword),
    CONSTRAINT fk_pdb_keyword_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);
//
DROP TABLE IF EXISTS pdb_chainid;
CREATE TABLE pdb_chainid(
    chain_id        VARCHAR(20) PRIMARY KEY,
    pdb_id          VARCHAR(10),
    model_no        INT,
    chain_no        VARCHAR(5),
    CONSTRAINT fk_chainid_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id)  ON DELETE CASCADE
);
//
CREATE TABLE pdb_chain(
    chain_id        VARCHAR(20),
    pdb_id          VARCHAR(10),
    model_no        INT,
    chain_no        VARCHAR(5),
    chain_seq       TEXT,
    CONSTRAINT pk_chain_pdb PRIMARY KEY (chain_id),
    CONSTRAINT fk_chain_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id)  ON DELETE CASCADE
);
//
DROP TABLE IF EXISTS meta_chain_pdb;
CREATE TABLE meta_chain_pdb(
    chain_id        VARCHAR(20),
    pdb_file        VARCHAR(60),
    CONSTRAINT pk_pdb_chain PRIMARY KEY (chain_id),
    CONSTRAINT fk_pdb_chain FOREIGN KEY (chain_id)
        REFERENCES pdb_chainid(chain_id)
);
//
DROP TABLE IF EXISTS pdb_chain_seq;
CREATE TABLE pdb_chain_seq(
    chain_id        VARCHAR(20),
    pro_len         INT,
    chain_seq       TEXT,
    CONSTRAINT pk_chain_seq PRIMARY KEY (chain_id)
);
//
DROP TABLE IF EXISTS meta_chain_faa;
CREATE TABLE meta_chain_faa(
    chain_id        VARCHAR(20),
    faa_file        VARCHAR(60) NOT NULL,
    UNIQUE (chain_id),
    CONSTRAINT fk_faa_chain FOREIGN KEY (chain_id)
        REFERENCES pdb_chainid(chain_id)
);
//
DROP TABLE IF EXISTS meta_combo_pdb;
CREATE TABLE meta_combo_pdb(
    pdb_id          VARCHAR(10),
    model_no        INT,
    chain_combo     VARCHAR(10),
    pdb_file        VARCHAR(60),
    CONSTRAINT fk_combo_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id)
);
//
DROP TABLE IF EXISTS meta_faa;
CREATE TABLE meta_faa(
    pdb_id          VARCHAR(10),
    faa_file        VARCHAR(60),
    UNIQUE (pdb_id),
    CONSTRAINT fk_faa_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id)
);
//
DROP TABLE IF EXISTS pdb_chain_type;
CREATE TABLE pdb_chain_type(
    pdb_id          VARCHAR(10),
    chain_id        VARCHAR(20),
    chain_type      VARCHAR(20),
    water_res       INT,
    total_res       INT,
    non_water_res   INT,
    aa_res          INT,
    dna_res         INT,
    rna_res         INT,
    ion_res         INT,
    hetero_res      INT,
    other_res       INT,
    UNIQUE (chain_id),
    CONSTRAINT fk_chain_type FOREIGN KEY (chain_id)
        REFERENCES pdb_chainid(chain_id)
);

"""
#####################################################################################

if __name__ == "__main__":
    bc = BuildComplex(verbose=True)
    for query in SQL_TEXT.split('//'): 
        res = bc.create_table(query)
        if res is None:
            print(f"ERROR. Check the above SQL: {query}\n\n")

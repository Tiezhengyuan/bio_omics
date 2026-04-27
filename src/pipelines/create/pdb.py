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

'''
pdb -> meta_mmcif
    -> meta_pdb
                -> pdb_compound
                -> pdb_method
                -> pdb_keyword
                -> pdb_chainid
                -> pdb_chain_seq
                -> pdb_faa
                -> meta_chain_faa
                -> meta_chain_pdb
                -> meta_combo2
                -> meta_combo3
'''

SQL_TEXT = """
CREATE TABLE pdb(
    pdb_id          VARCHAR(10),
    pdb_url         VARCHAR(50) NOT NULL,
    release_date    DATE NOT NULL,
    name            VARCHAR(500) NOT NULL,
    CONSTRAINT pk_pdb PRIMARY KEY (pdb_id)
);
//
DROP TABLE IF EXISTS meta_pdb;
CREATE TABLE meta_pdb(
    pdb_id          VARCHAR(10),
    relative_path   VARCHAR(60),
    pdb_group       INT,
    UNIQUE (pdb_id),
    CONSTRAINT fk_meta_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);
//
DROP TABLE IF EXISTS meta_mmcif;
CREATE TABLE meta_mmcif(
    pdb_id          VARCHAR(10),
    relative_path   VARCHAR(60) NOT NULL,
    pdb_group       INT,
    UNIQUE (pdb_id),
    CONSTRAINT fk_meta_mmcif FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);
//
DROP TABLE IF EXISTS pdb_compound;
CREATE TABLE pdb_compound(
    pdb_id          VARCHAR(10),
    compound_no     INT,
    chain_no        VARCHAR(5) COLLATE utf8mb4_bin NOT NULL,
    compound_text   TEXT,
    CONSTRAINT pk_compound PRIMARY KEY (pdb_id, compound_no, chain_no),
    CONSTRAINT fk_compound_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);
//
CREATE TABLE pdb_method(
    pdb_id              VARCHAR(10),
    structure_method    VARCHAR(100),
    resolution          FLOAT,
    CONSTRAINT pk_method PRIMARY KEY (pdb_id),
    CONSTRAINT fk_method_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
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
    chain_id        VARCHAR(20) COLLATE utf8mb4_bin NOT NULL,
    pdb_id          VARCHAR(10),
    model_no        INT,
    chain_no        VARCHAR(5) COLLATE utf8mb4_bin NOT NULL,
    idx             INT,
    UNIQUE(chain_id),
    CONSTRAINT fk_chainid_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

//
DROP TABLE IF EXISTS pdb_chain_seq;
CREATE TABLE pdb_chain_seq(
    pdb_id          VARCHAR(10),
    chain_id        VARCHAR(20) COLLATE utf8mb4_bin NOT NULL,
    first_chain_id  VARCHAR(20),
    pro_len         INT,
    chain_seq       TEXT,
    UNIQUE(chain_id),
    CONSTRAINT fk_chainseq_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);
//
DROP TABLE IF EXISTS meta_chain_faa;
CREATE TABLE meta_chain_faa(
    pdb_id          VARCHAR(10),
    chain_id        VARCHAR(20) COLLATE utf8mb4_bin NOT NULL,
    relative_path   VARCHAR(60),
    UNIQUE (pdb_id, chain_id),
    CONSTRAINT fk_faa_chain FOREIGN KEY (chain_id)
        REFERENCES pdb_chainid(chain_id) ON DELETE CASCADE
);
//
DROP TABLE IF EXISTS meta_chain_pdb;
CREATE TABLE meta_chain_pdb(
    pdb_id          VARCHAR(10),
    chain_id        VARCHAR(20) COLLATE utf8mb4_bin NOT NULL,
    chain_pdb       VARCHAR(60),
    UNIQUE (pdb_id, chain_id),
    CONSTRAINT fk_pdb_chain FOREIGN KEY (chain_id)
        REFERENCES pdb_chainid(chain_id) ON DELETE CASCADE
);
//
DROP TABLE IF EXISTS combo2;
CREATE TABLE combo2 (
    pdb_id          VARCHAR(10),
    combo_id        VARCHAR(40) COLLATE utf8mb4_bin NOT NULL,
    chain_combo     VARCHAR(10),
    chain_pdbs      VARCHAR(120),
    combo_pdb       VARCHAR(60),
    combo_csv       VARCHAR(60),
    combo_fasta     VARCHAR(60),
    UNIQUE(combo_id),
    CONSTRAINT fk_combo2_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);
//
DROP TABLE IF EXISTS combo3;
CREATE TABLE combo3 (
    pdb_id          VARCHAR(10),
    combo_id        VARCHAR(40) COLLATE utf8mb4_bin NOT NULL,
    chain_combo     VARCHAR(10),
    chain_pdbs      VARCHAR(120),
    combo_pdb       VARCHAR(60),
    combo_csv       VARCHAR(60),
    combo_fasta     VARCHAR(60),
    UNIQUE(combo_id),
    CONSTRAINT fk_combo3_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);
//
DROP TABLE IF EXISTS combo4;
CREATE TABLE combo4 (
    pdb_id          VARCHAR(10),
    combo_id        VARCHAR(40) COLLATE utf8mb4_bin NOT NULL,
    chain_combo     VARCHAR(10),
    chain_pdbs      VARCHAR(200),
    combo_pdb       VARCHAR(60),
    combo_csv       VARCHAR(60),
    combo_fasta     VARCHAR(60),
    UNIQUE(combo_id),
    CONSTRAINT fk_combo4_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);
//
DROP TABLE IF EXISTS combo5;
CREATE TABLE combo5 (
    pdb_id          VARCHAR(10),
    combo_id        VARCHAR(40) COLLATE utf8mb4_bin NOT NULL,
    chain_combo     VARCHAR(10),
    chain_pdbs      VARCHAR(300),
    combo_pdb       VARCHAR(60),
    combo_csv       VARCHAR(60),
    combo_fasta     VARCHAR(60),
    UNIQUE(combo_id),
    CONSTRAINT fk_combo5_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);
//
DROP TABLE IF EXISTS pdb_faa;
CREATE TABLE pdb_faa(
    pdb_id          VARCHAR(10),
    model_no        INT,
    relative_faa    VARCHAR(60),
    UNIQUE (pdb_id, model_no),
    CONSTRAINT fk_pdb_faa FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);
//
DROP TABLE IF EXISTS label_chain;
CREATE TABLE label_chain(
    pdb_id          VARCHAR(10),
    chain_id        VARCHAR(20) COLLATE utf8mb4_bin NOT NULL,
    first_chain_id  VARCHAR(20),
    chain_label     VARCHAR(20),
    category        VARCHAR(20),
    UNIQUE(chain_id, category),
    CONSTRAINT fk_label_chain_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);
//
DROP TABLE IF EXISTS label_complex;
CREATE TABLE label_complex(
    pdb_id          VARCHAR(10),
    model_no        INT,
    combo_id        VARCHAR(40) COLLATE utf8mb4_bin NOT NULL,
    chain_combo     VARCHAR(20),
    chain_label     VARCHAR(50),
    combo           INT,
    chain_ids       VARCHAR(100),
    category        VARCHAR(20),
    UNIQUE(combo_id, category),
    CONSTRAINT fk_label_complex_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);
//
DROP TABLE IF EXISTS pdb_protein;
CREATE TABLE pdb_protein (
    pdb_id          VARCHAR(10),
    relative_fasta  VARCHAR(60),
    pro_id          VARCHAR(20),
    chain           VARCHAR(40),
    auth_chain      VARCHAR(10),
    pro_len         INT,
    pro_seq         TEXT,
    UNIQUE(pdb_id, chain),
    CONSTRAINT fk_pdb_protein_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);
"""
#####################################################################################

if __name__ == "__main__":
    bc = BuildComplex(verbose=True)
    for query in SQL_TEXT.split('//'): 
        res = bc.create_table(query)
        if res is None:
            print(f"ERROR. Check the above SQL: {query}\n\n")



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

# deprecated
'''
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
'''

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
CREATE TABLE pdb_compound(
    pdb_id          VARCHAR(10),
    compound_no     INT,
    chains          VARCHAR(200),
    compound_text   TEXT,
    CONSTRAINT pk_compound PRIMARY KEY (pdb_id, compound_no),
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
    chain_id        VARCHAR(20),
    pdb_id          VARCHAR(10),
    model_no        INT,
    chain_no        VARCHAR(5),
    idx             INT,
    UNIQUE(chain_id),
    CONSTRAINT fk_chainid_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

//
DROP TABLE IF EXISTS pdb_chain_seq;
CREATE TABLE pdb_chain_seq(
    pdb_id          VARCHAR(10),
    chain_id        VARCHAR(20),
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
    chain_id        VARCHAR(20),
    relative_path   VARCHAR(60),
    UNIQUE (pdb_id, chain_id),
    CONSTRAINT fk_faa_chain FOREIGN KEY (chain_id)
        REFERENCES pdb_chainid(chain_id) ON DELETE CASCADE
);
//
DROP TABLE IF EXISTS meta_chain_pdb;
CREATE TABLE meta_chain_pdb(
    pdb_id          VARCHAR(10),
    chain_id        VARCHAR(20),
    chain_pdb       VARCHAR(60),
    UNIQUE (pdb_id, chain_id),
    CONSTRAINT fk_pdb_chain FOREIGN KEY (chain_id)
        REFERENCES pdb_chainid(chain_id) ON DELETE CASCADE
);
//
DROP TABLE IF EXISTS meta_combo2;
CREATE TABLE meta_combo2 (
    pdb_id          VARCHAR(10),
    combo_id        VARCHAR(40),
    chain_combo     VARCHAR(10),
    chain_pdbs      VARCHAR(120),
    combo_pdb       VARCHAR(60),
    combo_csv       VARCHAR(60),
    combo_fasta     VARCHAR(60),
    CONSTRAINT fk_combo2 FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);
//
DROP TABLE IF EXISTS meta_combo3;
CREATE TABLE meta_combo3 (
    pdb_id          VARCHAR(10),
    combo_id        VARCHAR(40),
    chain_combo     VARCHAR(10),
    chain_pdbs      VARCHAR(120),
    combo_pdb       VARCHAR(60),
    combo_csv       VARCHAR(60),
    combo_fasta     VARCHAR(60),
    CONSTRAINT fk_combo3 FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);
# //
# DROP TABLE IF EXISTS meta_combo2_pdb;
# CREATE TABLE meta_combo2_pdb(
#     pdb_id          VARCHAR(10),
#     model_no        INT,
#     combo_id        VARCHAR(40),
#     chain_combo     VARCHAR(10),
#     bounded_pdb     VARCHAR(60),
#     chain_pdbs      VARCHAR(60),
#     CONSTRAINT fk_combo2_pdb FOREIGN KEY (pdb_id)
#         REFERENCES pdb(pdb_id) ON DELETE CASCADE
# );
# //
# DROP TABLE IF EXISTS meta_combo2_seq;
# CREATE TABLE meta_combo2_seq(
#     pdb_id          VARCHAR(10),
#     model_no        INT,
#     combo_id        VARCHAR(40),
#     chain_combo     VARCHAR(10),
#     relative_csv    VARCHAR(60),
#     CONSTRAINT fk_combo2_seq FOREIGN KEY (pdb_id)
#         REFERENCES pdb(pdb_id) ON DELETE CASCADE
# );
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
DROP TABLE IF EXISTS pdb_complex;
CREATE TABLE pdb_complex(
    pdb_id          VARCHAR(10),
    combo_id        VARCHAR(40),
    chain_combo     VARCHAR(20),
    chain_label     VARCHAR(100),
    chains_num      INT,
    category        VARCHAR(20),
    CONSTRAINT fk_pdb_complex_pdb FOREIGN KEY (pdb_id)
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



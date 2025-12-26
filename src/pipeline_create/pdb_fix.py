#!/usr/bin/python
'''
requirements:
    - a empty database is created.
example:
    - python src/pipeline_create/pdb_fix.py
functions:
    - build tables
'''
from create_helper import *

'''
meta_pdb -> fixed_pdb -> fixed_chainid
                      -> fixed_chain
                      -> fixed_comb2
                      -> fixed_comb3
         -> meta_chain_pdb -> pdb_chain_fixed
'''

SQL_TEXT = """
DROP TABLE IF EXISTS fixed_pdb;
CREATE TABLE fixed_pdb (
    pdb_id          VARCHAR(10) NOT NULL,
    relative_path   VARCHAR(40),
    time0           DATETIME,
    created_at      TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    UNIQUE (pdb_id),
    CONSTRAINT fk_fixed_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id)
);
//
DROP TABLE IF EXISTS fixed_chainid;
CREATE TABLE fixed_chainid (
    pdb_id          VARCHAR(10),
    model_no        INT,
    fixed_model_no  INT,
    chain_no        VARCHAR(5),
    fixed_chain_no  VARCHAR(5),
    CONSTRAINT fk_chainid_fixed FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id)
);
//
DROP TABLE IF EXISTS fixed_chain;
CREATE TABLE fixed_chain (
    pdb_id          VARCHAR(10),
    chain_id        VARCHAR(20),
    model_no        INT,
    chain_no        VARCHAR(5),
    relative_pdb    VARCHAR(60),
    UNIQUE (pdb_id, chain_id),
    CONSTRAINT fk_fixed_chain FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id)
);
//
DROP TABLE IF EXISTS fixed_combo2;
CREATE TABLE fixed_combo2 (
    pdb_id          VARCHAR(10),
    model_no        INT,
    combo_id        VARCHAR(60),
    chain_combo     VARCHAR(10),
    relative_pdb    VARCHAR(60),
    UNIQUE (pdb_id, combo_id),
    CONSTRAINT fk_fixed_combo2 FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id)
);
//
DROP TABLE IF EXISTS fixed_combo3;
CREATE TABLE fixed_combo3 (
    pdb_id          VARCHAR(10),
    model_no        INT,
    combo_id        VARCHAR(60),
    chain_combo     VARCHAR(10),
    relative_pdb    VARCHAR(60),
    UNIQUE (pdb_id, combo_id),
    CONSTRAINT fk_fixed_combo3 FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id)
);

//
DROP TABLE IF EXISTS pdb_chain_fixed;
CREATE TABLE pdb_chain_fixed (
    chain_id        VARCHAR(20) PRIMARY KEY,
    relative_pdb    VARCHAR(60),
    time0           DATETIME,
    created_at      TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

"""
#####################################################################################

if __name__ == "__main__":
    bc = BuildComplex(verbose=True)
    for query in SQL_TEXT.split('//'): 
        res = bc.create_table(query)
        if res is None:
            print(f"ERROR. Check the above SQL: {query}\n\n")

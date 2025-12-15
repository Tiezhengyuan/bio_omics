#!/usr/bin/python
'''
requirements:
    - a empty database is created.
example:
    - python src/pipeline_create/pdb_simulation.py
functions:
    - build tables
'''
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)


from bioomics import BuildComplex

'''
relax chain/combo of pdb separately 
meta_pdb -> fixed_pdb -> fixed_chain -> relaxed_chain_amber14
                      -> fixed_comb2
                      -> fixed_comb3
'''

SQL_TEXT = """
DROP TABLE IF EXISTS fixed_pdb;
CREATE TABLE fixed_pdb (
    pdb_id          VARCHAR(10),
    pdb_file        VARCHAR(40),
    created_at      TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    UNIQUE (pdb_id),
    CONSTRAINT fk_fixed_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id)
);
//
DROP TABLE IF EXISTS fixed_chain;
CREATE TABLE fixed_chain (
    chain_id        VARCHAR(20) PRIMARY KEY,
    pdb_id          VARCHAR(10),
    model_no        INT,
    chain_no        VARCHAR(5),
    pdb_file        VARCHAR(60),
    CONSTRAINT fk_fixed_chain FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id)
);
//
DROP TABLE IF EXISTS relaxed_chain_amber14;
CREATE TABLE relaxed_chain_amber14 (
    chain_id            VARCHAR(20),
    original_energy     VARCHAR(50),
    relaxed_energy      VARCHAR(50),
    endpoint            VARCHAR(80),
    created_at          TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    UNIQUE (chain_id),
    CONSTRAINT fk_relaxed_chain_amber14 FOREIGN KEY (chain_id)
        REFERENCES fixed_chain(chain_id)
);
//
DROP TABLE IF EXISTS fixed_combo2;
CREATE TABLE fixed_combo2 (
    complex_id      INT PRIMARY KEY AUTO_INCREMENT,
    pdb_id          VARCHAR(10),
    model_no        INT,
    chain_combo     VARCHAR(10),
    pdb_file        VARCHAR(60),
    UNIQUE (pdb_id, model_no, chain_combo),
    CONSTRAINT fk_fixed_combo2 FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id)
);
//
DROP TABLE IF EXISTS fixed_combo3;
CREATE TABLE fixed_combo3 (
    complex_id      INT PRIMARY KEY AUTO_INCREMENT,
    pdb_id          VARCHAR(10),
    model_no        INT,
    chain_combo     VARCHAR(10),
    pdb_file        VARCHAR(60),
    UNIQUE (pdb_id, model_no, chain_combo),
    CONSTRAINT fk_fixed_combo3 FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id)
);

"""
#####################################################################################

if __name__ == "__main__":
    bc = BuildComplex(verbose=True)
    for query in SQL_TEXT.split('//'): 
        res = bc.create_table(query)
        if res is None:
            print(f"ERROR. Check the above SQL: {query}\n\n")

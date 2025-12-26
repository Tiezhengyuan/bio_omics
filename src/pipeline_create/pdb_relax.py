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
# method 1:
fixed_pdb -> relaxed_pdb_amber99sb
          -> relaxed_pdb_amber14 -> relaxed_amber14_chain
                                 -> relaxed_amber14_combo2

# method 2: relax chain/combo of pdb separately 
fixed_pdb -> fixed_chain -> relaxed_chain_amber14

# method 3: chain/combo -> fixed + relaxed
meta_chain_pdb -> meta_chain_relaxed
meta_combo2_pdb -> meta_combo2_relaxed

'''

SQL_TEXT = """
DROP TABLE IF EXISTS relaxed_pdb_amber99sb;
CREATE TABLE relaxed_pdb_amber99sb (
    pdb_id              VARCHAR(10),
    original_energy     VARCHAR(50),
    relaxed_energy      VARCHAR(50),
    relative_pdb        VARCHAR(50),
    time0               DATETIME,
    created_at          TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    UNIQUE (pdb_id),
    CONSTRAINT fk_relaxed_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id)
);
//
DROP TABLE IF EXISTS relaxed_pdb_amber14;
CREATE TABLE relaxed_pdb_amber14 (
    pdb_id              VARCHAR(10),
    original_energy     VARCHAR(50),
    relaxed_energy      VARCHAR(50),
    relative_pdb        VARCHAR(50),
    time0               DATETIME,
    created_at          TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    UNIQUE (pdb_id),
    CONSTRAINT fk_relaxed_pdb_amber14 FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id)
);
//
DROP TABLE IF EXISTS relaxed_amber14_chain;
CREATE TABLE relaxed_amber14_chain (
    chain_id            VARCHAR(20),
    pdb_id              VARCHAR(10),
    model_no            INT,
    chain_no            VARCHAR(5),
    relative_path       VARCHAR(50),
    UNIQUE (chain_id),
    CONSTRAINT fk_relaxed_amber14_chain FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id)
);
//
DROP TABLE IF EXISTS relaxed_amber14_combo2;
CREATE TABLE relaxed_amber14_combo2 (
    pdb_id              VARCHAR(10),
    combo_id            VARCHAR(40),
    model_no            INT,
    chain_combo         VARCHAR(20),
    endpoint_prefix     VARCHAR(80),
    UNIQUE (combo_id),
    CONSTRAINT fk_relaxed_amber14_combo2 FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id)
);
//
DROP TABLE IF EXISTS relaxed_chain_amber14;
CREATE TABLE relaxed_chain_amber14 (
    chain_id            VARCHAR(20),
    original_energy     VARCHAR(50),
    relaxed_energy      VARCHAR(50),
    endpoint            VARCHAR(80),
    time0               DATETIME,
    created_at          TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    UNIQUE (chain_id)
);
"""
#####################################################################################

if __name__ == "__main__":
    bc = BuildComplex(verbose=True)
    for query in SQL_TEXT.split('//'): 
        res = bc.create_table(query)
        if res is None:
            print(f"ERROR. Check the above SQL: {query}\n\n")

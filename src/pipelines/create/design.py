#!/usr/bin/python
'''
requirements:
    - a empty database is created.
example:
    - python src/pipeline_create/design.py
functions:
    - build tables
'''
from create_helper import *



SQL_TEXT = """
DROP TABLE IF EXISTS design_dimer;
CREATE TABLE design_dimer (
    design_id       VARCHAR(100) COLLATE utf8mb4_bin PRIMARY KEY,
    category        VARCHAR(10),
    chain_labels    VARCHAR(20),
    target_name     VARCHAR(20),
    target_id       VARCHAR(40),
    binder_id       VARCHAR(40),
    af_faa          VARCHAR(100) NOT NULL,
    target_pdb      VARCHAR(100),
    binder_pdb      VARCHAR(100)
);
//
DROP TABLE IF EXISTS design_multimer;
CREATE TABLE design_multimer (
    design_id       VARCHAR(100) COLLATE utf8mb4_bin PRIMARY KEY,
    category        VARCHAR(10),
    chain_labels    VARCHAR(20),
    combo           INT,
    target_name     VARCHAR(20),
    target_id       VARCHAR(40),
    binder_id       VARCHAR(100),
    af_faa          VARCHAR(100) NOT NULL,
    target_pdb      VARCHAR(100),
    binder_pdb      VARCHAR(200)
);
//
DROP TABLE IF EXISTS predict;
CREATE TABLE predict (
    predictor      VARCHAR(10),
    design_id      VARCHAR(100),
    host           VARCHAR(10),
    endpoint       VARCHAR(100),
    start_at       TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    end_at         DATETIME,
    status         VARCHAR(10),
    UNIQUE(predictor, design_id)
);
"""
#####################################################################################

if __name__ == "__main__":
    bc = BuildComplex(verbose=True)
    for query in SQL_TEXT.split('//'): 
        res = bc.create_table(query)
        if res is None:
            print(f"ERROR. Check the above SQL: {query}\n\n")

#!/usr/bin/python
'''
requirements:
    - a empty database is created.
example:
    - python src/pipeline_create/predict.py
functions:
    - build tables
'''
from create_helper import *

SQL_TEXT = """
DROP TABLE IF EXISTS chain_alphafold;
CREATE TABLE chain_alphafold (
    chain_id            VARCHAR(20),
    endpoint            VARCHAR(60),
    predicted_at        DATEIME,
    created_at          TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    UNIQUE (chain_id)
);
//
DROP TABLE IF EXISTS combo2_alphafold;
CREATE TABLE combo2_alphafold (
    combo_id            VARCHAR(20),
    endpoint            VARCHAR(60),
    predicted_at        DATEIME,
    created_at          TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    UNIQUE (combo_id)
);

"""
#####################################################################################

if __name__ == "__main__":
    bc = BuildComplex(verbose=True)
    for query in SQL_TEXT.split('//'): 
        res = bc.create_table(query)
        if res is None:
            print(f"ERROR. Check the above SQL: {query}\n\n")

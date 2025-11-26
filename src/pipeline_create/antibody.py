#!/usr/bin/python
'''
requirements:
    - an empty database is created.
example: python src/pipelines/create_ab.py
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
DROP TABLE IF EXISTS ab_aacdb;
CREATE TABLE ab_aacdb (
    pdb_id          VARCHAR(10),
    antibody_chains VARCHAR(20),
    antigen_chains  VARCHAR(20),
    antibody        VARCHAR(100),
    protein	        VARCHAR(200),
    drug	        VARCHAR(500),
    targets	        VARCHAR(200),
    method          VARCHAR(50),
    organism        VARCHAR(200),
    resolution      FLOAT,
    reference       VARCHAR(500)
);
//
DROP TABLE IF EXISTS meta_sabdab;
CREATE TABLE meta_sabdab (
    pdb_id          VARCHAR(10),
    summary_file    VARCHAR(200),
    UNIQUE(pdb_id)
);
//
DROP TABLE IF EXISTS ab_sabdab;
CREATE TABLE ab_sabdab (
    pdb_id          VARCHAR(100),
    model           INT,
    heavy_chain     VARCHAR(20),
    light_chain     VARCHAR(20),
    antigen_chain   VARCHAR(30),
    antigen_type    VARCHAR(50),
    resolution      FLOAT
);
//
DROP TABLE IF EXISTS ab_abybank;
CREATE TABLE ab_abybank (
    pdb_id          VARCHAR(10),
    chain_type      VARCHAR(20),
    chain_numbering VARCHAR(20)
);
"""

#####################################################################################

if __name__ == "__main__":
    bc = BuildComplex(verbose=True)
    for query in SQL_TEXT.split('//'): 
        res = bc.create_table(query)
        if res is None:
            print(f"ERROR. Check the above SQL: {query}\n\n")

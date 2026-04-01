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
DROP TABLE IF EXISTS meta_execute;
CREATE TABLE meta_execute(
    table_name      VARCHAR(50),
    pdb_group       INT,
    host            VARCHAR(10),
    ip_addr         VARCHAR(30),
    script          VARCHAR(100),
    start_at        TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    end_at          DATETIME,
    UNIQUE (table_name, pdb_group, host)
);
//
DROP TABLE IF EXISTS meta;
CREATE TABLE meta(
    pdb_group       INT,
    host            VARCHAR(10),
    ip_addr         VARCHAR(30),
    script          VARCHAR(100),
    start_at        TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    end_at          DATETIME,
    UNIQUE (table_name, pdb_group, host)
);
"""
#####################################################################################

if __name__ == "__main__":
    bc = BuildComplex(verbose=True)
    for query in SQL_TEXT.split('//'): 
        res = bc.create_table(query)
        if res is None:
            print(f"ERROR. Check the above SQL: {query}\n\n")



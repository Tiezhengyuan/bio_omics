#!/usr/bin/python
'''
requirements:
    - an empty database is created.
example: python src/pipelines/create_protein.py
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
DROP TABLE IF EXISTS pro_pep;
CREATE TABLE pro_pep (
    pep_id      INT PRIMARY KEY,
    seq         VARCHAR(10)
);
//
DROP TABLE IF EXISTS protein_seq;
CREATE TABLE protein_seq (
    pro_id      INT PRIMARY KEY,
    prefix      VARCHAR(3),
    pro_len     INT,
    pro_seq     TEXT
);
CREATE INDEX index_prolen on protein_seq (pro_len);
//
DROP TABLE IF EXISTS protein_source;
CREATE TABLE protein_source (
    pro_id          INT,
    source_name     VARCHAR(50),
    source_id       VARCHAR(50),
    UNIQUE(source_name, source_id),
    CONSTRAINT fk_protein_source FOREIGN KEY (pro_id)
        REFERENCES protein_seq(pro_id) ON DELETE CASCADE
);
"""

#####################################################################################

if __name__ == "__main__":
    bc = BuildComplex(verbose=True)
    for query in SQL_TEXT.split('//'): 
        res = bc.create_table(query)
        if res is None:
            print(f"ERROR. Check the above SQL: {query}\n\n")

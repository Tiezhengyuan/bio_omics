#!/usr/bin/python
'''
requirements:
    - a empty database is created.
example: python src/pipelines/create_imgt.py
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
DROP TABLE IF EXISTS meta_bfactor;
CREATE TABLE meta_bfactor(
    pdb_id          VARCHAR(20),
    pickle_file     VARCHAR(100),
    CONSTRAINT pk_meta_bfactor PRIMARY KEY (pdb_id),
    CONSTRAINT fk_bfactor_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);
//
DROP TABLE IF EXISTS pdb_bfactor;
CREATE TABLE pdb_bfactor(
    pdb_id          VARCHAR(20),
    model_no        INT,
    chain_no        VARCHAR(5),
    avg_bfactor     FLOAT
);
//
DROP TABLE IF EXISTS pdb_secondary;
CREATE TABLE pdb_secondary(
    pdb_id          VARCHAR(10),
    model_no        INT,
    chain_no        VARCHAR(5),
    dssp_seq        TEXT,
    CONSTRAINT fk_secondary_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id)
);
//
DROP TABLE IF EXISTS meta_pdb_distance;
CREATE TABLE meta_pdb_distance(
    pdb_id          VARCHAR(20),
    pickle_file     VARCHAR(100),
    UNIQUE (pdb_id),
    CONSTRAINT fk_distance_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id)
);
//
DROP TABLE IF EXISTS meta_pdb_sasa;
CREATE TABLE meta_pdb_sasa(
    pdb_id          VARCHAR(20),
    pickle_file     VARCHAR(100),
    UNIQUE (pdb_id),
    CONSTRAINT fk_sasa_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id)
);
//
DROP TABLE IF EXISTS meta_pdb_igblastp;
CREATE TABLE meta_pdb_igblastp(
    pdb_id          VARCHAR(20),
    aln_file        VARCHAR(100) NOT NULL,
    UNIQUE (pdb_id),
    CONSTRAINT fk_igblastp_pdb FOREIGN KEY (pdb_id)
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

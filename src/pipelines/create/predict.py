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

'''
meta_pdb -> view_pdb_chain -> chain_alphafold -> chain_alphafold_plddt
                                              -> chain_model_measure
         -> meta_combo2 -> combo2_alphafold

'''


SQL_TEXT = """
DROP TABLE IF EXISTS chain_alphafold;
CREATE TABLE chain_alphafold (
    chain_id            VARCHAR(20),
    endpoint            VARCHAR(60),
    predicted_at        DATETIME,
    created_at          TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    UNIQUE (chain_id)
);
DROP TABLE IF EXISTS chain_plddt;
CREATE TABLE chain_plddt (
    predict_id      INT AUTO_INCREMENT PRIMARY KEY,
    chain_id        VARCHAR(20),
    ranking         INT,
    predictor       VARCHAR(10),
    avg_plddt       FLOAT,
    relative_pkl    VARCHAR(100),
    relative_pdb    VARCHAR(100),
    UNIQUE (chain_id, relative_pdb)
);
//
DROP TABLE IF EXISTS chain_rmsd;
CREATE TABLE chain_rmsd (
    pdb_id          VARCHAR(10),
    chain_id        VARCHAR(20),
    rmsd            FLOAT,
    chain_status    VARCHAR(20),
    native_pdb      VARCHAR(100),
    model_pdb       VARCHAR(100),
    UNIQUE(native_pdb, model_pdb),
    CONSTRAINT fk_chain_rmsd_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);
DROP TABLE IF EXISTS chain_tmalign;
CREATE TABLE chain_tmalign(
    pdb_id          VARCHAR(10),
    chain_id        VARCHAR(20),
    tm1             FLOAT,
    tm2             FLOAT,
    rmsd            FLOAT,
    chain_status    VARCHAR(20),
    native_pdb      VARCHAR(100),
    model_pdb       VARCHAR(100),
    UNIQUE(native_pdb, model_pdb),
    CONSTRAINT fk_chain_tmalign_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);
//
DROP TABLE IF EXISTS combo2_alphafold;
CREATE TABLE combo2_alphafold (
    combo_id            VARCHAR(20),
    endpoint            VARCHAR(60),
    predicted_at        DATETIME,
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

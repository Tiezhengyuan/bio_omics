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
meta_pdb -> view_pdb_chain -> chain_af -> chain_af_plddt
                                              -> chain_model_measure
         -> meta_combo2 -> combo2_af
         -> meta_combo3 -> combo3_af

'''


SQL_TEXT = """
DROP TABLE IF EXISTS chain_af;
CREATE TABLE chain_af (
    pdb_id              VARCHAR(10),
    chain_id            VARCHAR(20) COLLATE utf8mb4_bin NOT NULL,
    endpoint            VARCHAR(60),
    host                VARCHAR(10),
    predicted_at        DATETIME,
    created_at          TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    UNIQUE (chain_id),
    CONSTRAINT fk_chain_alpha_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);
DROP TABLE IF EXISTS chain_af_summary;
CREATE TABLE chain_af_summary (
    chain_id        VARCHAR(20) COLLATE utf8mb4_bin NOT NULL,
    ranking         INT,
    avg_plddt       FLOAT,
    avg_ptm         FLOAT,
    max_pae         FLOAT,
    relative_json   VARCHAR(150),
    relaxed_pdb     VARCHAR(150),
    unrelaxed_pdb   VARCHAR(150),
    UNIQUE (chain_id, relaxed_pdb)
);
//
DROP TABLE IF EXISTS chain_af_plddt;
CREATE TABLE chain_af_plddt (
    chain_id        VARCHAR(20) COLLATE utf8mb4_bin NOT NULL,
    ranking         INT,
    avg_plddt       FLOAT,
    relative_pkl    VARCHAR(200),
    relative_pdb    VARCHAR(100),
    UNIQUE (chain_id, relative_pdb)
);
//
DROP TABLE IF EXISTS chain_rmsd;
CREATE TABLE chain_rmsd (
    pdb_id          VARCHAR(10),
    chain_id        VARCHAR(20) COLLATE utf8mb4_bin NOT NULL,
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
    chain_id        VARCHAR(20) COLLATE utf8mb4_bin NOT NULL,
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
DROP TABLE IF EXISTS combo2_af;
CREATE TABLE combo2_af (
    pdb_id              VARCHAR(10),
    combo_id            VARCHAR(20) COLLATE utf8mb4_bin NOT NULL,
    endpoint            VARCHAR(60),
    host                VARCHAR(10),
    predicted_at        DATETIME,
    created_at          TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    UNIQUE (combo_id),
    CONSTRAINT fk_combo2_alpha_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);
DROP TABLE IF EXISTS combo3_af;
CREATE TABLE combo3_af (
    pdb_id              VARCHAR(10),
    combo_id            VARCHAR(20) COLLATE utf8mb4_bin NOT NULL,
    endpoint            VARCHAR(60),
    host                VARCHAR(10),
    predicted_at        DATETIME,
    created_at          TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    UNIQUE (combo_id),
    CONSTRAINT fk_combo3_alpha_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);
"""
#####################################################################################

if __name__ == "__main__":
    bc = BuildComplex(verbose=True)
    for query in SQL_TEXT.split('//'): 
        res = bc.create_table(query)
        if res is None:
            print(f"ERROR. Check the above SQL: {query}\n\n")

#!/usr/bin/python
'''
requirements:
    - a empty database is created.
example: python src/pipelines/create_uniprot.py
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
CREATE TABLE uniprot (
    uniprot_acc     VARCHAR(20),
    name            VARCHAR(100),
    organism        VARCHAR(200),
    pro_len         INT,
    pro_seq         TEXT,
    CONSTRAINT pk_uniprot PRIMARY KEY (uniprot_acc)
);
//
CREATE TABLE uniprot_keyword (
    uniprot_acc     VARCHAR(20),
    keyword         VARCHAR(300),
    UNIQUE (uniprot_acc, keyword),
    CONSTRAINT fk_uniprot_keyword_uniprot FOREIGN KEY (uniprot_acc)
        REFERENCES uniprot(uniprot_acc) ON DELETE CASCADE
);
//
CREATE TABLE uniprot_dbxrefs (
    uniprot_acc     VARCHAR(20),
    db_name         VARCHAR(50),
    db_acc          VARCHAR(50),
    UNIQUE (uniprot_acc, db_name, db_acc),
    CONSTRAINT fk_uniprot_dbxrefs_uniprot FOREIGN KEY (uniprot_acc)
        REFERENCES uniprot(uniprot_acc) ON DELETE CASCADE
);
//
DROP TABLE IF EXISTS uniprot_feature;
CREATE TABLE uniprot_feature (
    uniprot_acc     VARCHAR(20),
    feature         VARCHAR(100),
    start           VARCHAR(50),
    end             VARCHAR(50),
    seq             TEXT,
    note            TEXT,
    evidence        TEXT,
    CONSTRAINT fk_feature_uniprot FOREIGN KEY (uniprot_acc)
        REFERENCES uniprot(uniprot_acc) ON DELETE CASCADE
);
//
DROP TABLE IF EXISTS uniprot_pro_seq;
CREATE TABLE uniprot_pro_seq (
    uniprot_acc  VARCHAR(20) PRIMARY KEY,
    proseq_name  VARCHAR(20),
    seq_id       INT
);
//
DROP TABLE IF EXISTS uniprot_af;
CREATE TABLE uniprot_af (
    uniprot_acc     VARCHAR(20),
    pdb_path        VARCHAR(200),
    CONSTRAINT fk_alphafold_uniprot FOREIGN KEY (uniprot_acc)
    REFERENCES uniprot(uniprot_acc)
);
//
DROP TABLE IF EXISTS uniprot_af_plddt;
CREATE TABLE uniprot_af_plddt (
    uniprot_acc     VARCHAR(20),
    model_no        INT,
    chain_no        VARCHAR(5),
    avg_plddt       FLOAT,
    pickle_file     VARCHAR(200),
    CONSTRAINT fk_plddt_uniprot FOREIGN KEY (uniprot_acc)
    REFERENCES uniprot(uniprot_acc)
);
//
DROP TABLE IF EXISTS uniprot_af_rmsd;
CREATE TABLE uniprot_af_rmsd (
    chain_id        VARCHAR(20),
    uniprot_acc     VARCHAR(20),
    ca_rmsd         FLOAT,
    native_pdb      VARCHAR(100),
    model_pdb       VARCHAR(100),
    UNIQUE (chain_id, uniprot_acc)
);
"""

#####################################################################################

if __name__ == "__main__":
    bc = BuildComplex(verbose=True)
    for query in SQL_TEXT.split('//'): 
        res = bc.create_table(query)
        if res is None:
            print(f"ERROR. Check the above SQL: {query}\n\n")

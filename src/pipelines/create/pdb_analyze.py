#!/usr/bin/python
'''
requirements:
    - a empty database is created.
example: python src/pipeline_create/pdb_analyze.py
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

'''
meta_pdb -> meta_bfactor -> pdb_bfactor
         -> pdb_chainid -> pdb_chain_type
                        -> pdb_chain_pdb -> pdb_dssp
         -> meta_combo2 -> combo2_distance
                        -> combo2_sasa
        -> pdb_distance

'''

SQL_TEXT = """
DROP TABLE IF EXISTS meta_bfactor;
CREATE TABLE meta_bfactor(
    pdb_id           VARCHAR(10),
    relative_pkl     VARCHAR(60),
    CONSTRAINT pk_meta_bfactor PRIMARY KEY (pdb_id),
    CONSTRAINT fk_bfactor_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);
//
DROP TABLE IF EXISTS pdb_bfactor;
CREATE TABLE pdb_bfactor(
    pdb_id          VARCHAR(10),
    model_no        INT,
    chain_no        VARCHAR(5) COLLATE utf8mb4_bin NOT NULL,
    avg_bfactor     FLOAT,
    bfactor         TEXT,
    UNIQUE(pdb_id, model_no, chain_no),
    CONSTRAINT fk_pdb_bfactor_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);
//
DROP TABLE IF EXISTS pdb_chain_type;
CREATE TABLE pdb_chain_type(
    pdb_id          VARCHAR(10),
    chain_id        VARCHAR(20),
    chain_type      VARCHAR(20),
    water_res       INT,
    total_res       INT,
    non_water_res   INT,
    aa_res          INT,
    dna_res         INT,
    rna_res         INT,
    ion_res         INT,
    hetero_res      INT,
    other_res       INT,
    UNIQUE (chain_id),
    CONSTRAINT fk_chain_type_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);
//
DROP TABLE IF EXISTS chain_dssp;
CREATE TABLE chain_dssp(
    pdb_id    VARCHAR(10),
    chain_id  VARCHAR(20) COLLATE utf8mb4_bin NOT NULL,
    model_no  INT,
    chain_no  VARCHAR(5),
    dssp_seq  TEXT,
    host      VARCHAR(10),
    UNIQUE(chain_id),
    CONSTRAINT fk_chaindssp_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);
//
DROP TABLE IF EXISTS combo2_distance;
CREATE TABLE combo2_distance(
    pdb_id          VARCHAR(10),
    combo_id        VARCHAR(40),
    relative_pkl    VARCHAR(60),
    min_dist        FLOAT,
    min_res1        VARCHAR(10),
    min_res2        VARCHAR(10),
    CONSTRAINT fk_combo2_distance_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);
//
DROP TABLE IF EXISTS combo2_sasa;
CREATE TABLE combo2_sasa(
    pdb_id           VARCHAR(10),
    combo_id         VARCHAR(40),
    delta_total_sasa FLOAT,
    relative_pkl     VARCHAR(60),
    CONSTRAINT fk_combo2_sasa_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);
//
DROP TABLE IF EXISTS combo_contacts;
CREATE TABLE combo_contacts(
    pdb_id                  VARCHAR(10),
    combo_id                VARCHAR(40),
    binding_affinity        FLOAT,
    dissociation_constant   FLOAT,
    relative_prodigy        VARCHAR(100),
    CONSTRAINT fk_combo_contacts_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);
//
DROP TABLE IF EXISTS equal_seq_rmsd;
CREATE TABLE equal_seq_rmsd(
    chain_id        VARCHAR(10),
    equal_chain_id  VARCHAR(10),
    rmsd            FLOAT,
    chain_pdb       VARCHAR(60),
    equal_chain_pdb VARCHAR(60),
    UNIQUE(chain_id, equal_chain_id)
);
//
DROP TABLE IF EXISTS dist_amber14;
CREATE TABLE dist_amber14(
    pdb_id          VARCHAR(10),
    combo_id        VARCHAR(40),
    relative_pkl    VARCHAR(60),
    CONSTRAINT fk_dist_amber14_pdb FOREIGN KEY (pdb_id)
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

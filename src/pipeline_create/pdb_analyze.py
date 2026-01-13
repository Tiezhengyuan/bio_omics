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
         -> meta_combo2_pdb -> meta_pdb_distance
                            -> meta_pdb_sasa
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
    chain_no        VARCHAR(5),
    avg_bfactor     FLOAT,
    bfactor         TEXT,
    CONSTRAINT fk_chain_bfactor_pdb FOREIGN KEY (pdb_id)
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
DROP TABLE IF EXISTS pdb_dssp;
CREATE TABLE pdb_dssp(
    pdb_id          VARCHAR(10),
    model_no        INT,
    chain_no        VARCHAR(5),
    dssp_seq        TEXT,
    CONSTRAINT fk_secondary_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);
//
DROP TABLE IF EXISTS meta_pdb_distance;
CREATE TABLE meta_pdb_distance(
    pdb_id          VARCHAR(20),
    combo_id        VARCHAR(40),
    relative_pkl    VARCHAR(60),
    CONSTRAINT fk_meta_distance_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);
//
DROP TABLE IF EXISTS meta_pdb_sasa;
CREATE TABLE meta_pdb_sasa(
    pdb_id          VARCHAR(20),
    combo_id        VARCHAR(40),
    relative_pkl    VARCHAR(60),
    CONSTRAINT fk_sasa_pdb FOREIGN KEY (pdb_id)
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

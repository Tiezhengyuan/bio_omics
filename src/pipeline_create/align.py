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

'''
align_imgt: igblastp of pdb chains on IMGT V-region
align_vregion: FR and CDR idetity in v-region

pdb_faa -> meta_blastp5  -> align_uniprot
        -> meta_igblastp -> align_imgt
                         -> align_vregion

'''

SQL_TEXT = """
DROP TABLE IF EXISTS meta_blastp5;
CREATE TABLE meta_blastp5(
    pdb_id          VARCHAR(10),
    model_no        INT,
    relative_aln    VARCHAR(60),
    CONSTRAINT fk_blastp_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id)
);
//
DROP TABLE IF EXISTS align_uniprot;
CREATE TABLE align_uniprot(
    pdb_id          VARCHAR(10),
    chain_id        VARCHAR(20),
    uniprot_acc     VARCHAR(20),
    evalue          DOUBLE,
    bit_score       INT,
    identity        FLOAT,
    UNIQUE (chain_id)
);
//
DROP TABLE IF EXISTS meta_igblastp;
CREATE TABLE meta_igblastp(
    pdb_id          VARCHAR(10),
    model_no        INT,
    relative_aln    VARCHAR(60),
    CONSTRAINT fk_meta_igblastp_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id)
);
//
DROP TABLE IF EXISTS align_vregion;
CREATE TABLE align_vregion (
    pdb_id          VARCHAR(10),
    chain_id        VARCHAR(20),
    allele_name     VARCHAR(20),
    specie          VARCHAR(40),
    gene_name       VARCHAR(20),
    gene_family     VARCHAR(20),
    chain_type      VARCHAR(10),
    isotype         VARCHAR(10),
    evalue          DOUBLE,
    bit_score       INT,
    region_names    VARCHAR(100),
    UNIQUE(chain_id),
    CONSTRAINT fk_imgt_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id)
);
//
DROP TABLE IF EXISTS align_vfrag;
CREATE TABLE align_vfrag (
    pdb_id          VARCHAR(10),
    chain_id        VARCHAR(20),
    region_name     VARCHAR(20),
    identity        FLOAT,
    seq_from        INT,
    seq_to          INT,
    seq             TEXT,
    UNIQUE(chain_id, region_name),
    CONSTRAINT fk_vregion_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id)
);
//
DROP TABLE IF EXISTS imgt_inn_blastp;
CREATE TABLE imgt_inn_blastp(
    query_pro_id    INT,
    inn_chain_id    VARCHAR(20),
    ref_pro_id      INT,
    pdb_chain_id    VARCHAR(20),
    evalue          DOUBLE,
    bit_score       INT,
    identity        FLOAT
);
"""

#####################################################################################

if __name__ == "__main__":
    bc = BuildComplex(verbose=True)
    for query in SQL_TEXT.split('//'): 
        res = bc.create_table(query)
        if res is None:
            print(f"ERROR. Check the above SQL: {query}\n\n")

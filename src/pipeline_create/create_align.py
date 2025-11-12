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
'''

SQL_TEXT = """
DROP TABLE IF EXISTS meta_blastp5;
CREATE TABLE meta_blastp5(
    pdb_id          VARCHAR(20),
    aln_file        VARCHAR(100) NOT NULL,
    CONSTRAINT fk_blastp_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id)
);
//
DROP TABLE IF EXISTS align_uniprot;
CREATE TABLE align_uniprot(
    chain_id        VARCHAR(20),
    pdb_id          VARCHAR(10),
    uniprot_acc     VARCHAR(20),
    evalue          DOUBLE,
    bit_score       INT,
    identity        FLOAT,
    CONSTRAINT pk_align_uniprot PRIMARY KEY (chain_id)
);
//
DROP TABLE IF EXISTS align_imgt;
CREATE TABLE align_imgt(
    chain_id        VARCHAR(20),
    allele_name     VARCHAR(20),
    gene_name       VARCHAR(20),
    gene_family     VARCHAR(20),
    chain_type      VARCHAR(10),
    isotype         VARCHAR(10),
    evalue          DOUBLE,
    bit_score       INT,
    CONSTRAINT pk_imgt PRIMARY KEY (chain_id),
    CONSTRAINT fk_imgt_chain FOREIGN KEY (chain_id)
        REFERENCES pdb_chain(chain_id)
);
//
DROP TABLE IF EXISTS align_vregion;
CREATE TABLE align_vregion (
    chain_id        VARCHAR(20),
    pdb_id          VARCHAR(10),
    region_name     VARCHAR(20),
    identity        FLOAT,
    start           INT,
    end             INT,
    seq             TEXT,
    CONSTRAINT pk_vregion PRIMARY KEY (chain_id, region_name),
    CONSTRAINT fk_vregion_chain FOREIGN KEY (chain_id)
        REFERENCES pdb_chain(chain_id)
);
"""

#####################################################################################

if __name__ == "__main__":
    bc = BuildComplex(verbose=True)
    for query in SQL_TEXT.split('//'): 
        res = bc.create_table(query)
        if res is None:
            print(f"ERROR. Check the above SQL: {query}\n\n")

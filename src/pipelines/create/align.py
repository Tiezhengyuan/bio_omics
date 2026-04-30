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
        -> meta_vdj -> align_imgt
                    -> align_vregion
                    -> align_dregion
                    -> align_jregion

'''

SQL_TEXT = """
DROP TABLE IF EXISTS meta_blastp5;
CREATE TABLE meta_blastp5(
    pdb_id          VARCHAR(10),
    model_no        INT,
    relative_aln    VARCHAR(60),
    CONSTRAINT fk_blastp_pdb FOREIGN KEY (pdb_id) REFERENCES pdb(pdb_id)
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
DROP TABLE IF EXISTS meta_vdj;
CREATE TABLE meta_vdj(
    pdb_id      VARCHAR(10),
    model_no    INT,
    v_aln       VARCHAR(60),
    d_aln       VARCHAR(60),
    j_aln       VARCHAR(60),
    CONSTRAINT fk_align_vdj_pdb FOREIGN KEY (pdb_id) REFERENCES pdb(pdb_id)
);
//
DROP TABLE IF EXISTS align_vfrag;
CREATE TABLE align_vfrag (
    pdb_id          VARCHAR(10),
    chain_id        VARCHAR(20),
    fragment        VARCHAR(20),
    identity        FLOAT,
    seq_from        INT,
    seq_to          INT,
    seq             TEXT,
    UNIQUE(chain_id, fragment),
    CONSTRAINT fk_vfrag_pdb FOREIGN KEY (pdb_id) REFERENCES pdb(pdb_id)
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
    CONSTRAINT fk_imgt_pdb FOREIGN KEY (pdb_id) REFERENCES pdb(pdb_id)
);
//
DROP TABLE IF EXISTS align_dregion;
CREATE TABLE align_dregion (
    pdb_id          VARCHAR(10),
    chain_id        VARCHAR(20),
    query_len       INT,
    hit_id          VARCHAR(50),
    organism        VARCHAR(30),
    evalue          DOUBLE,
    bit_score       INT,
    match_seq       VARCHAR(50),
    identities      INT,
    align_len       INT,
    query_seq       VARCHAR(50),
    query_start     INT,
    query_end       INT,
    sbjct_seq       VARCHAR(50),
    sbjct_start     INT,
    sbjct_end       INT,
    CONSTRAINT fk_dregion_pdb FOREIGN KEY (pdb_id) REFERENCES pdb(pdb_id)
);
//
DROP TABLE IF EXISTS align_jregion;
CREATE TABLE align_jregion (
    pdb_id          VARCHAR(10),
    chain_id        VARCHAR(20),
    query_len       INT,
    allele_name     VARCHAR(20),
    specie          VARCHAR(40),
    gene_name       VARCHAR(20),
    gene_family     VARCHAR(20),
    chain_type      VARCHAR(10),
    isotype         VARCHAR(10),
    evalue          DOUBLE,
    bit_score       INT,
    match_seq       VARCHAR(50),
    identities      INT,
    align_len       INT,
    query_seq       VARCHAR(50),
    query_start     INT,
    query_end       INT,
    sbjct_seq       VARCHAR(50),
    sbjct_start     INT,
    sbjct_end       INT,
    CONSTRAINT fk_jregion_pdb FOREIGN KEY (pdb_id) REFERENCES pdb(pdb_id)
);


"""

#####################################################################################

if __name__ == "__main__":
    bc = BuildComplex(verbose=True)
    for query in SQL_TEXT.split('//'): 
        res = bc.create_table(query)
        if res is None:
            print(f"ERROR. Check the above SQL: {query}\n\n")

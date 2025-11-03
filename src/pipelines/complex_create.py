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

# PDB data
pdb = """
CREATE TABLE pdb(
    pdb_id          VARCHAR(10),
    pdb_url         VARCHAR(50) NOT NULL,
    release_date    DATE NOT NULL,
    name            VARCHAR(500) NOT NULL,
    CONSTRAINT pk_pdb PRIMARY KEY (pdb_id)
);
//
CREATE TABLE pdb_method(
    pdb_id              VARCHAR(10),
    structure_method    VARCHAR(100),
    resolution          FLOAT,
    CONSTRAINT pk_method PRIMARY KEY (pdb_id),
    CONSTRAINT fk_method_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id)  ON DELETE CASCADE
);
//
CREATE TABLE pdb_compound(
    pdb_id          VARCHAR(10),
    compound_no     INT,
    chains          VARCHAR(200),
    compound_text   TEXT,
    CONSTRAINT pk_compound PRIMARY KEY (pdb_id, compound_no),
    CONSTRAINT fk_compound_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id)  ON DELETE CASCADE
);
//
CREATE TABLE pdb_keyword(
    pdb_id          VARCHAR(10),
    keyword         VARCHAR(200),
    UNIQUE (pdb_id, keyword),
    CONSTRAINT fk_pdb_keyword_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);
//
CREATE TABLE pdb_chain(
    chain_id        VARCHAR(20),
    pdb_id          VARCHAR(10),
    model_no        INT,
    chain_no        VARCHAR(5),
    chain_seq       TEXT,
    proseq_id       INT,
    CONSTRAINT pk_chain_pdb PRIMARY KEY (chain_id),
    CONSTRAINT fk_chain_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id)  ON DELETE CASCADE
);
//
DROP TABLE IF EXISTS pdb_proseq;
CREATE TABLE pdb_proseq(
    proseq_id       INT PRIMARY KEY AUTO_INCREMENT,
    chunk_no        INT,
    pro_seq         TEXT
);
//
CREATE TABLE pdb_bfactor(
    pdb_id          VARCHAR(10),
    model_no        INT,
    chain_no        VARCHAR(5),
    avg_bfactor     FLOAT,
    CONSTRAINT fk_bfactor_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id)
);
//
CREATE TABLE pdb_secondary(
    pdb_id          VARCHAR(10),
    model_no        INT,
    chain_no        VARCHAR(5),
    dssp_seq        TEXT,
    CONSTRAINT fk_secondary_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id)
);

"""

#####################
meta = """
DROP TABLE IF EXISTS meta_pdb
CREATE TABLE meta_pdb(
    pdb_id          VARCHAR(10),
    pdb_file        VARCHAR(100) NOT NULL,
    pdb_group       INT,
    UNIQUE (pdb_id),
    CONSTRAINT fk_meta_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id)  ON DELETE CASCADE
);
//
DROP TABLE IF EXISTS meta_chain_faa;
CREATE TABLE meta_chain_faa(
    chain_id        VARCHAR(20),
    faa_file        VARCHAR(100) NOT NULL,
    UNIQUE (chain_id),
    CONSTRAINT fk_faa_chain FOREIGN KEY (chain_id)
        REFERENCES pdb_chain(chain_id)
);
//
DROP TABLE IF EXISTS meta_chain_pdb;
CREATE TABLE meta_chain_pdb(
    pdb_id          VARCHAR(20),
    model_no        INT,
    chain_combo     VARCHAR(10),
    pdb_file        VARCHAR(100) NOT NULL,
    CONSTRAINT fk_pdb_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id)
);
//
DROP TABLE IF EXISTS meta_faa;
CREATE TABLE meta_faa(
    pdb_id          VARCHAR(20),
    faa_file        VARCHAR(100) NOT NULL,
    UNIQUE (pdb_id),
    CONSTRAINT fk_faa_pdb FOREIGN KEY (pdb_id)
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


##########################
# uniprot data: accession primary key
uniprot = """
CREATE TABLE uniprot (
    uniprot_acc     VARCHAR(20),
    name            VARCHAR(100),
    organism        VARCHAR(200),
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
DROP TABLE IF EXISTS uniprot_meta_alphafold;
CREATE TABLE uniprot_meta_alphafold (
    uniprot_acc     VARCHAR(20),
    pdb_path        VARCHAR(200),
    CONSTRAINT fk_alphafold_uniprot FOREIGN KEY (uniprot_acc)
    REFERENCES uniprot(uniprot_acc)
);
//
DROP TABLE IF EXISTS uniprot_alphafold_plddt;
CREATE TABLE uniprot_alphafold_plddt (
    uniprot_acc     VARCHAR(20),
    model_no        INT,
    chain_no        VARCHAR(5),
    avg_plddt       FLOAT,
    pickle_file     VARCHAR(200),
    CONSTRAINT fk_plddt_uniprot FOREIGN KEY (uniprot_acc)
    REFERENCES uniprot(uniprot_acc)
);
//
DROP TABLE IF EXISTS uniprot_rmsd;
CREATE TABLE uniprot_rmsd (
    chain_id        VARCHAR(20),
    native_pdb      VARCHAR(200),
    uniprot_acc     VARCHAR(20),
    model_pdb       VARCHAR(200),
    ca_rmsd         FLOAT,
    UNIQUE (chain_id, uniprot_acc)
);

"""

'''
align_imgt: igblastp of pdb chains on IMGT V-region
align_vregion: FR and CDR idetity in v-region
'''
antibody = """
CREATE TABLE imgt_genelist (
    species             VARCHAR(50),
    gene_db             VARCHAR(20),
    gene_function       VARCHAR(20),
    gene_definition     VARCHAR(200),
    number_alleles      VARCHAR(5),
    chromosome          VARCHAR(5),
    chromosomal_loc     VARCHAR(20),
    ligm_db_ref_gene    VARCHAR(100),
    ncbi_gene           VARCHAR(20),
    chain               VARCHAR(10),
    domain              VARCHAR(5)
);
//
DROP TABLE IF EXISTS imgt_gene;
CREATE TABLE imgt_gene (
    germline_id         INT PRIMARY KEY,
    specie              VARCHAR(100),
    isotype             VARCHAR(10),
    chain_type          VARCHAR(10),
    gene_name           VARCHAR(20),
    allele_name         VARCHAR(20),
    region_name         VARCHAR(50)
);
//
DROP TABLE IF EXISTS imgt_geneseq;
CREATE TABLE imgt_geneseq (
    germline_id     INT,
    seq_type        VARCHAR(100),
    seq             TEXT,
    UNIQUE (germline_id, seq_type),
    CONSTRAINT fk_imgt_genseq_gene FOREIGN KEY (germline_id)
    REFERENCES imgt_gene(germline_id) ON DELETE CASCADE
);
//

"""

align = """
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

# views
view = """
DROP VIEW IF EXISTS view_pdb;
CREATE VIEW view_pdb AS
    SELECT P.pdb_id, P.pdb_url, B.pdb_file, P.release_date,
        P.name, M.structure_method, M.resolution
    FROM pdb P, meta_pdb B, pdb_method M
    WHERE P.pdb_id = B.pdb_id
        AND P.pdb_id = M.pdb_id
        AND M.resolution IS NOT NULL
;
//
DROP VIEW IF EXISTS view_pdb_chain;
CREATE VIEW view_pdb_chain AS
    SELECT C.chain_id, C.pdb_id, C.model_no, C.chain_no,
        M.structure_method, M.resolution, B.avg_bfactor,
        C.chain_seq, F.faa_file, P.pdb_file
    FROM pdb_chain C
    LEFT JOIN meta_chain_faa F
        ON C.chain_id = F.chain_id
    LEFT JOIN pdb_method M
        ON C.pdb_id = M.pdb_id
    LEFT JOIN pdb_bfactor B
        ON C.pdb_id = B.pdb_id AND C.model_no = B.model_no
            AND C.chain_no = B.chain_no
    LEFT JOIN meta_chain_pdb P
        ON C.pdb_id = P.pdb_id AND C.model_no = P.model_no
            AND C.chain_no = P.chain_combo
;
//
DROP VIEW IF EXISTS view_antibody;
CREATE VIEW view_antibody AS
    SELECT I.allele_name, I.gene_name, I.gene_family, I.chain_type,
        C.pdb_id, C.chain_id, C.model_no, C.chain_no, C.chain_seq,
        C.structure_method, C.resolution, C.faa_file, C.pdb_file
    FROM view_pdb_chain C, align_imgt I
    WHERE C.chain_id = I.chain_id
        AND I.evalue <= 1e-10
        AND I.bit_score >= 100
        AND I.isotype = 'IG'
;
//



DROP VIEW IF EXISTS view_pdb_uniprot;
CREATE VIEW view_pdb_uniprot AS
    SELECT A.chain_id, A.pdb_id, A.uniprot_acc, U.name,
        K.keyword, C.chain_no, C.pro_seq
    FROM align_uniprot A
    LEFT JOIN uniprot U
        ON A.uniprot_acc = U.uniprot_acc
    LEFT JOIN uniprot_keyword K
        ON A.uniprot_acc = K.uniprot_acc
    LEFT JOIN pdb_chain C
        ON A.chain_id = C.chain_id
    WHERE A.bit_score > 100
        AND A.evalue < 1e-10
        AND C.model_no = 0
        AND C.pro_seq IS NOT NULL;
//

DROP VIEW IF EXISTS view_antibody_complex;
CREATE VIEW view_antibody_complex AS
    SELECT * FROM view_pdb_uniprot
    WHERE pdb_id IN(
        SELECT DISTINCT pdb_id
        FROM view_pdb_uniprot
        WHERE keyword = 'Immunoglobulin'
            OR keyword LIKE '%antibody%'
    );

"""


#####################################################################################

if __name__ == "__main__":
    db_name = 'COMPLEX2'
    bc = BuildComplex(verbose=True)
    for sql in (pdb, meta, uniprot, align, antibody, view):
        for query in sql.split('//'): 
            res = bc.create_table(query)
            if res is None:
                print(f"ERROR. Check the above SQL: {query}\n\n")

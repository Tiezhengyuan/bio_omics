/*
align_imgt: igblastp of pdb chains on IMGT V-region
align_vregion: FR and CDR idetity in v-region

pdb_faa -> meta_blastp5  -> align_uniprot
        -> meta_vdj -> align_imgt
                    -> align_vregion -> label_chain
                                     -> label_complex
                    -> align_dregion
                    -> align_jregion
*/

DROP TABLE IF EXISTS meta_blastp5;
CREATE TABLE meta_blastp5(
    pdb_id          VARCHAR(10),
    model_no        INT,
    relative_aln    VARCHAR(60),
    CONSTRAINT fk_blastp_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

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

DROP TABLE IF EXISTS align_vdj;
CREATE TABLE align_vdj(
    pdb_id      VARCHAR(10),
    model_no    INT,
    v_aln       VARCHAR(60),
    d_aln       VARCHAR(60),
    j_aln       VARCHAR(60),
    CONSTRAINT fk_align_vdj_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

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
    CONSTRAINT fk_vfrag_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

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
    CONSTRAINT fk_vregion_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id)  ON DELETE CASCADE
);

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
    CONSTRAINT fk_dregion_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

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
    CONSTRAINT fk_jregion_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS label_chain;
CREATE TABLE label_chain(
    pdb_id          VARCHAR(10),
    chain_id        VARCHAR(20) COLLATE utf8mb4_bin NOT NULL,
    first_chain_id  VARCHAR(20),
    chain_label     VARCHAR(20),
    category        VARCHAR(20),
    UNIQUE(chain_id, category),
    CONSTRAINT fk_label_chain_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS label_complex;
CREATE TABLE label_complex(
    pdb_id          VARCHAR(10),
    model_no        INT,
    combo_id        VARCHAR(40) COLLATE utf8mb4_bin NOT NULL,
    chain_combo     VARCHAR(20),
    chain_label     VARCHAR(50),
    combo           INT,
    chain_ids       VARCHAR(100),
    category        VARCHAR(20),
    UNIQUE(combo_id, category),
    CONSTRAINT fk_label_complex_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);
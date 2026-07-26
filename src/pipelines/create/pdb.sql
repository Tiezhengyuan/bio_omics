/*
pdb -> meta_mmcif
    -> meta_pdb
                -> pdb_compound
                -> pdb_method
                -> pdb_keyword
                -> pdb_protein
                -> meta_bfactor -> pdb_bfactor

                -> pdb_chainid
                -> pdb_chain_seq + proseq_*-> pdb_pro_seq
                -> pdb_faa
                -> chain_faa
                -> pdb_chain_type
                -> chain_pdb -> chain_dssp

*/


/* pdb */
CREATE TABLE pdb(
    pdb_id          VARCHAR(10),
    pdb_url         VARCHAR(50) NOT NULL,
    release_date    DATE NOT NULL,
    name            VARCHAR(500) NOT NULL,
    CONSTRAINT pk_pdb PRIMARY KEY (pdb_id)
);

DROP TABLE IF EXISTS meta_pdb;
CREATE TABLE meta_pdb(
    pdb_id          VARCHAR(10),
    relative_pdb    VARCHAR(60),
    relative_cif    VARCHAR(60),
    pdb_group       INT,
    UNIQUE (pdb_id)
);

DROP TABLE IF EXISTS pdb_compound;
CREATE TABLE pdb_compound(
    pdb_id          VARCHAR(10),
    compound_no     INT,
    chain_no        VARCHAR(5) COLLATE utf8mb4_bin NOT NULL,
    compound_text   TEXT,
    CONSTRAINT pk_compound PRIMARY KEY (pdb_id, compound_no, chain_no),
    CONSTRAINT fk_compound_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

CREATE TABLE pdb_method(
    pdb_id              VARCHAR(10),
    structure_method    VARCHAR(100),
    resolution          FLOAT,
    CONSTRAINT pk_method PRIMARY KEY (pdb_id),
    CONSTRAINT fk_method_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

CREATE TABLE pdb_keyword(
    pdb_id          VARCHAR(10),
    keyword         VARCHAR(200),
    UNIQUE (pdb_id, keyword),
    CONSTRAINT fk_pdb_keyword_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS pdb_protein;
CREATE TABLE pdb_protein (
    pdb_id          VARCHAR(10),
    relative_fasta  VARCHAR(60),
    pro_id          VARCHAR(20),
    chain           VARCHAR(80) COLLATE utf8mb4_bin,
    auth_chain      VARCHAR(10),
    pro_len         INT GENERATED ALWAYS AS (LENGTH(pro_seq)) STORED,
    pro_seq         TEXT,
    UNIQUE(pdb_id, chain),
    CONSTRAINT fk_pdb_protein_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS meta_bfactor;
CREATE TABLE meta_bfactor(
    pdb_id           VARCHAR(10),
    relative_pkl     VARCHAR(60),
    CONSTRAINT pk_meta_bfactor PRIMARY KEY (pdb_id),
    CONSTRAINT fk_bfactor_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

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


/* pdb_chain */
DROP TABLE IF EXISTS pdb_chainid;
CREATE TABLE pdb_chainid(
    chain_id        VARCHAR(20) COLLATE utf8mb4_bin NOT NULL,
    pdb_id          VARCHAR(10),
    model_no        INT,
    chain_no        VARCHAR(5) COLLATE utf8mb4_bin NOT NULL,
    idx             INT,
    UNIQUE(chain_id),
    CONSTRAINT fk_chainid_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS pdb_chain_seq;
CREATE TABLE pdb_chain_seq(
    pdb_id          VARCHAR(10),
    chain_id        VARCHAR(20) COLLATE utf8mb4_bin NOT NULL,
    first_chain_id  VARCHAR(20),
    pro_len         INT GENERATED ALWAYS AS (LENGTH(chain_seq)) STORED,
    chain_seq       TEXT,
    UNIQUE(chain_id),
    CONSTRAINT fk_chainseq_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS pdb_faa;
CREATE TABLE pdb_faa(
    pdb_id          VARCHAR(10),
    model_no        INT,
    relative_faa    VARCHAR(60),
    UNIQUE (pdb_id, model_no),
    CONSTRAINT fk_pdb_faa FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS chain_faa;
CREATE TABLE chain_faa(
    pdb_id          VARCHAR(10),
    chain_id        VARCHAR(20) COLLATE utf8mb4_bin NOT NULL,
    relative_path   VARCHAR(60),
    UNIQUE (pdb_id, chain_id),
    CONSTRAINT fk_faa_chain FOREIGN KEY (chain_id)
        REFERENCES pdb_chainid(chain_id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS pdb_pro_seq;
CREATE TABLE pdb_pro_seq(
    pdb_id      VARCHAR(10),
    chain_id    VARCHAR(20) COLLATE utf8mb4_bin NOT NULL,
    proseq_name VARCHAR(20),
    seq_id      INT,
    UNIQUE (chain_id),
    CONSTRAINT fk_pdbproseq_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

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

DROP TABLE IF EXISTS chain_pdb;
CREATE TABLE chain_pdb(
    pdb_id          VARCHAR(10),
    chain_id        VARCHAR(20) COLLATE utf8mb4_bin NOT NULL,
    chain_pdb       VARCHAR(60),
    UNIQUE (chain_id),
    CONSTRAINT fk_pdb_chain FOREIGN KEY (chain_id)
        REFERENCES pdb_chainid(chain_id) ON DELETE CASCADE
);

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
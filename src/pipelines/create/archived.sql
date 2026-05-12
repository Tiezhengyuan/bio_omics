/*
method 1:
fixed_pdb -> relaxed_pdb_amber99sb
          -> relaxed_pdb_amber14 -> relaxed_amber14_chain
                                 -> relaxed_amber14_combo2
                                 -> relaxed_amber14_combo3

method 2: relax chain/combo of pdb separately 
fixed_pdb -> fixed_chain -> chain_relaxed_amber14

method 3: chain/combo -> fixed + relaxed

meta_pdb -> fixed_pdb -> fixed_chainid
                      -> fixed_chain
                      -> fixed_comb2
                      -> fixed_comb3
         -> meta_chain_pdb -> pdb_chain_fixed

*/

DROP TABLE IF EXISTS meta_mmcif;
CREATE TABLE meta_mmcif(
    pdb_id          VARCHAR(10),
    relative_path   VARCHAR(60) NOT NULL,
    pdb_group       INT,
    UNIQUE (pdb_id),
    CONSTRAINT fk_meta_mmcif FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS relaxed_pdb_amber99sb;
CREATE TABLE relaxed_pdb_amber99sb (
    pdb_id              VARCHAR(10),
    original_energy     VARCHAR(50),
    relaxed_energy      VARCHAR(50),
    relative_pdb        VARCHAR(50),
    time0               DATETIME,
    created_at          TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    UNIQUE (pdb_id),
    CONSTRAINT fk_relaxed_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS relaxed_pdb_amber14;
CREATE TABLE relaxed_pdb_amber14 (
    pdb_id              VARCHAR(10),
    original_energy     VARCHAR(50),
    relaxed_energy      VARCHAR(50),
    relative_pdb        VARCHAR(50),
    time0               DATETIME,
    created_at          TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    UNIQUE (pdb_id),
    CONSTRAINT fk_relaxed_pdb_amber14 FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS chain_relaxed_amber14;
CREATE TABLE chain_relaxed_amber14 (
    chain_id            VARCHAR(20) COLLATE utf8mb4_bin NOT NULL,
    original_energy     VARCHAR(50),
    relaxed_energy      VARCHAR(50),
    endpoint            VARCHAR(80),
    time0               DATETIME,
    created_at          TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    UNIQUE (chain_id)
);

DROP TABLE IF EXISTS relaxed_amber14_chain;
CREATE TABLE relaxed_amber14_chain (
    pdb_id              VARCHAR(10),
    chain_id            VARCHAR(20),
    model_no            INT,
    chain_no            VARCHAR(5),
    relative_pdb        VARCHAR(50),
    UNIQUE (pdb_id, chain_id),
    CONSTRAINT fk_relaxed_amber14_chain FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS relaxed_amber14_combo2;
CREATE TABLE relaxed_amber14_combo2 (
    pdb_id              VARCHAR(10),
    combo_id            VARCHAR(40),
    model_no            INT,
    chain_combo         VARCHAR(20),
    endpoint_prefix     VARCHAR(80),
    UNIQUE (combo_id),
    CONSTRAINT fk_relaxed_amber14_combo2 FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS relaxed_amber14_combo3;
CREATE TABLE relaxed_amber14_combo3 (
    pdb_id              VARCHAR(10),
    combo_id            VARCHAR(40),
    model_no            INT,
    chain_combo         VARCHAR(20),
    endpoint_prefix     VARCHAR(80),
    UNIQUE (combo_id),
    CONSTRAINT fk_relaxed_amber14_combo3 FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS fixed_pdb;
CREATE TABLE fixed_pdb (
    pdb_id          VARCHAR(10) NOT NULL,
    relative_path   VARCHAR(40),
    time0           DATETIME,
    created_at      TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    UNIQUE (pdb_id),
    CONSTRAINT fk_fixed_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS fixed_chainid;
CREATE TABLE fixed_chainid (
    pdb_id          VARCHAR(10),
    model_no        INT,
    fixed_model_no  INT,
    chain_no        VARCHAR(5),
    fixed_chain_no  VARCHAR(5),
    CONSTRAINT fk_chainid_fixed FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS fixed_chain;
CREATE TABLE fixed_chain (
    pdb_id          VARCHAR(10),
    chain_id        VARCHAR(20),
    model_no        INT,
    chain_no        VARCHAR(5),
    relative_pdb    VARCHAR(60),
    UNIQUE (pdb_id, chain_id),
    CONSTRAINT fk_fixed_chain FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS fixed_combo2;
CREATE TABLE fixed_combo2 (
    pdb_id          VARCHAR(10),
    model_no        INT,
    combo_id        VARCHAR(60),
    chain_combo     VARCHAR(10),
    relative_pdb    VARCHAR(60),
    UNIQUE (pdb_id, combo_id),
    CONSTRAINT fk_fixed_combo2 FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS fixed_combo3;
CREATE TABLE fixed_combo3 (
    pdb_id          VARCHAR(10),
    model_no        INT,
    combo_id        VARCHAR(60),
    chain_combo     VARCHAR(10),
    relative_pdb    VARCHAR(60),
    UNIQUE (pdb_id, combo_id),
    CONSTRAINT fk_fixed_combo3 FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);


DROP TABLE IF EXISTS pdb_chain_fixed;
CREATE TABLE pdb_chain_fixed (
    chain_id        VARCHAR(20) PRIMARY KEY,
    relative_pdb    VARCHAR(60),
    time0           DATETIME,
    created_at      TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

DROP TABLE IF EXISTS dist_amber14;
CREATE TABLE dist_amber14(
    pdb_id          VARCHAR(10),
    combo_id        VARCHAR(40),
    relative_pkl    VARCHAR(60),
    CONSTRAINT fk_dist_amber14_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS pro_pep;
CREATE TABLE pro_pep (
    pep_id      INT PRIMARY KEY,
    seq         VARCHAR(10)
);

DROP TABLE IF EXISTS protein_seq;
CREATE TABLE protein_seq (
    pro_id      INT PRIMARY KEY,
    prefix      VARCHAR(3),
    pro_len     INT,
    pro_seq     TEXT
);
CREATE INDEX index_prolen on protein_seq (pro_len);

DROP TABLE IF EXISTS protein_source;
CREATE TABLE protein_source (
    pro_id          INT,
    source_name     VARCHAR(50),
    source_id       VARCHAR(50),
    UNIQUE(source_name, source_id),
    CONSTRAINT fk_protein_source FOREIGN KEY (pro_id)
        REFERENCES protein_seq(pro_id) ON DELETE CASCADE
);
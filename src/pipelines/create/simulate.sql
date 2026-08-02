/*
meta_pdb -> pdb_chain_seq -> chain_pdb -> chain_relaxed
                                       -> equal_chain_rmsd
         -> combo2 -> combo2_relaxed
                   -> combo2_distance
         -> combo3 -> combo3_relaxed
         -> combo4 -> combo4_relaxed
         -> combo5 -> combo5_relaxed

                    ->combo_sasa
                    ->combo_contacts
*/

DROP TABLE IF EXISTS combo2;
CREATE TABLE combo2 (
    pdb_id          VARCHAR(10),
    combo_id        VARCHAR(40) COLLATE utf8mb4_bin NOT NULL,
    chain_combo     VARCHAR(10) NOT NULL,
    chain_pdbs      VARCHAR(120) NOT NULL,
    combo_pdb       VARCHAR(60) NOT NULL,
    combo_csv       VARCHAR(60) NOT NULL,
    combo_fasta     VARCHAR(60) NOT NULL,
    UNIQUE(combo_id),
    CONSTRAINT fk_combo2_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS combo3;
CREATE TABLE combo3 (
    pdb_id          VARCHAR(10),
    combo_id        VARCHAR(40) COLLATE utf8mb4_bin NOT NULL,
    chain_combo     VARCHAR(10) NOT NULL,
    chain_pdbs      VARCHAR(120) NOT NULL,
    combo_pdb       VARCHAR(60) NOT NULL,
    combo_csv       VARCHAR(60) NOT NULL,
    combo_fasta     VARCHAR(60) NOT NULL,
    UNIQUE(combo_id),
    CONSTRAINT fk_combo3_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS combo4;
CREATE TABLE combo4 (
    pdb_id          VARCHAR(10),
    combo_id        VARCHAR(40) COLLATE utf8mb4_bin NOT NULL,
    chain_combo     VARCHAR(10) NOT NULL,
    chain_pdbs      VARCHAR(200) NOT NULL,
    combo_pdb       VARCHAR(60) NOT NULL,
    combo_csv       VARCHAR(60) NOT NULL,
    combo_fasta     VARCHAR(60) NOT NULL,
    UNIQUE(combo_id),
    CONSTRAINT fk_combo4_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS combo5;
CREATE TABLE combo5 (
    pdb_id          VARCHAR(10),
    combo_id        VARCHAR(40) COLLATE utf8mb4_bin NOT NULL,
    chain_combo     VARCHAR(10) NOT NULL,
    chain_pdbs      VARCHAR(300) NOT NULL,
    combo_pdb       VARCHAR(60) NOT NULL,
    combo_csv       VARCHAR(60) NOT NULL,
    combo_fasta     VARCHAR(60) NOT NULL,
    UNIQUE(combo_id),
    CONSTRAINT fk_combo5_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS chain_fixed;
CREATE TABLE chain_fixed (
    pdb_id              VARCHAR(10),
    chain_id            VARCHAR(20) COLLATE utf8mb4_bin NOT NULL,
    fixed_pdb           VARCHAR(80) NOT NULL,
    UNIQUE (chain_id),
    CONSTRAINT fk_chain_fixed_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS chain_relaxed;
CREATE TABLE chain_relaxed (
    pdb_id              VARCHAR(10),
    chain_id            VARCHAR(20) COLLATE utf8mb4_bin NOT NULL,
    time0               DATETIME,
    relaxed_pdb         VARCHAR(80),
    original_energy     VARCHAR(60),
    relaxed_energy      VARCHAR(60),
    created_at          TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    UNIQUE (chain_id),
    CONSTRAINT fk_chain_relaxed_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS equal_seq_rmsd;
CREATE TABLE equal_seq_rmsd(
    chain_id        VARCHAR(10),
    equal_chain_id  VARCHAR(10),
    rmsd            FLOAT NOT NULL,
    chain_pdb       VARCHAR(60),
    equal_chain_pdb VARCHAR(60),
    UNIQUE(chain_id, equal_chain_id)
);

DROP TABLE IF EXISTS combo2_fixed;
CREATE TABLE combo2_fixed(
    pdb_id              VARCHAR(10),
    combo_id            VARCHAR(40) COLLATE utf8mb4_bin NOT NULL,
    fixed_pdb           VARCHAR(80) NOT NULL,
    UNIQUE (combo_id),
    CONSTRAINT fk_combo2_fixed_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS combo2_relaxed;
CREATE TABLE combo2_relaxed (
    pdb_id              VARCHAR(10),
    combo_id            VARCHAR(40) COLLATE utf8mb4_bin NOT NULL,
    time0               DATETIME,
    relaxed_pdb         VARCHAR(80),
    original_energy     VARCHAR(60),
    relaxed_energy      VARCHAR(60),
    created_at          TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    UNIQUE (combo_id),
    CONSTRAINT fk_combo2_relaxed_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS combo3_fixed;
CREATE TABLE combo3_fixed(
    pdb_id              VARCHAR(10),
    combo_id            VARCHAR(40) COLLATE utf8mb4_bin NOT NULL,
    fixed_pdb           VARCHAR(80) NOT NULL,
    UNIQUE (combo_id),
    CONSTRAINT fk_combo3_fixed_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS combo3_relaxed;
CREATE TABLE combo3_relaxed (
    pdb_id              VARCHAR(10),
    combo_id            VARCHAR(40) COLLATE utf8mb4_bin NOT NULL,
    time0               DATETIME,
    relaxed_pdb         VARCHAR(80),
    original_energy     VARCHAR(60),
    relaxed_energy      VARCHAR(60),
    created_at          TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    UNIQUE (combo_id),
    CONSTRAINT fk_combo3_relaxed_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS combo4_fixed;
CREATE TABLE combo4_fixed(
    pdb_id              VARCHAR(10),
    combo_id            VARCHAR(40) COLLATE utf8mb4_bin NOT NULL,
    fixed_pdb           VARCHAR(80) NOT NULL,
    UNIQUE (combo_id),
    CONSTRAINT fk_combo4_fixed_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS combo4_relaxed;
CREATE TABLE combo4_relaxed (
    pdb_id              VARCHAR(10),
    combo_id            VARCHAR(40) COLLATE utf8mb4_bin NOT NULL,
    time0               DATETIME,
    relaxed_pdb         VARCHAR(80),
    original_energy     VARCHAR(60),
    relaxed_energy      VARCHAR(60),
    created_at          TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    UNIQUE (combo_id),
    CONSTRAINT fk_combo4_relaxed_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS combo5_fixed;
CREATE TABLE combo5_fixed(
    pdb_id              VARCHAR(10),
    combo_id            VARCHAR(40) COLLATE utf8mb4_bin NOT NULL,
    fixed_pdb           VARCHAR(80) NOT NULL,
    UNIQUE (combo_id),
    CONSTRAINT fk_combo5_fixed_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS combo5_relaxed;
CREATE TABLE combo5_relaxed (
    pdb_id              VARCHAR(10),
    combo_id            VARCHAR(40) COLLATE utf8mb4_bin NOT NULL,
    time0               DATETIME,
    relaxed_pdb         VARCHAR(80),
    original_energy     VARCHAR(60),
    relaxed_energy      VARCHAR(60),
    created_at          TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    UNIQUE (combo_id),
    CONSTRAINT fk_combo5_relaxed_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS combo2_distance;
CREATE TABLE combo2_distance(
    pdb_id          VARCHAR(10),
    combo_id        VARCHAR(40),
    relative_pkl    VARCHAR(60),
    min_dist        FLOAT,
    min_res1        VARCHAR(10),
    min_res2        VARCHAR(10),
    created_at      TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    CONSTRAINT fk_combo2_distance_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS combo_sasa;
CREATE TABLE combo_sasa(
    pdb_id            VARCHAR(10),
    combo_id          VARCHAR(40),
    delta_total_sasa  FLOAT,
    relative_pkl      VARCHAR(60),
    created_at        TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    CONSTRAINT fk_combo_sasa_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS combo_contacts;
CREATE TABLE combo_contacts(
    pdb_id                  VARCHAR(10),
    combo_id                VARCHAR(40),
    binding_affinity        FLOAT,
    dissociation_constant   FLOAT,
    relative_prodigy        VARCHAR(100),
    created_at              TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    CONSTRAINT fk_combo_contacts_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);
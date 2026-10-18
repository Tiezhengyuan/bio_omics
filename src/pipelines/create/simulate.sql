/*
meta_pdb -> pdb_chain_seq -> chain_pdb -> equal_chain_rmsd
         -> combo2 -> combo2_distance
         -> combo3
         -> combo4
         -> combo5

                    ->combo_sasa
                    ->combo_contacts

chain_fixed, combo2..5_fixed, chain_relaxed, combo2..5_relaxed,
pdb_qc and pdb_relaxed_qc moved to structure_profiling/sql/simulate.sql
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

/* chain_fixed, chain_relaxed moved to structure_profiling/sql/simulate.sql */

DROP TABLE IF EXISTS equal_seq_rmsd;
CREATE TABLE equal_seq_rmsd(
    chain_id        VARCHAR(10),
    equal_chain_id  VARCHAR(10),
    rmsd            FLOAT NOT NULL,
    chain_pdb       VARCHAR(60),
    equal_chain_pdb VARCHAR(60),
    UNIQUE(chain_id, equal_chain_id)
);

/* combo2_fixed, combo2_relaxed moved to structure_profiling/sql/simulate.sql */

/* combo3_fixed, combo3_relaxed moved to structure_profiling/sql/simulate.sql */

/* combo4_fixed, combo4_relaxed moved to structure_profiling/sql/simulate.sql */

/*
combo5_fixed, combo5_relaxed, pdb_qc and pdb_relaxed_qc
moved to structure_profiling/sql/simulate.sql
*/
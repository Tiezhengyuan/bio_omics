/*
label_chain -> chain_af -> chain_af_summary -> chain_af_plddt
                                            -> chain_af_rmsd
                                            -> chain_tmalign
label_complex -> combo2 -> combo2_af
              -> combo3 -> combo3_af
              -> combo4 -> combo4_af
              -> combo5 -> combo5_af
*/

DROP TABLE IF EXISTS pdb_pred_params;
CREATE TABLE pdb_pred_params(
    pdb_id            VARCHAR(10),
    input_id          VARCHAR(20) COLLATE utf8mb4_bin NOT NULL,
    input_table_name  VARCHAR(20),
    input_col_name    VARCHAR(20),
    category          VARCHAR(20),
    predictor         VARCHAR(20),
    query_faa         VARCHAR(200) NOT NULL,
    outdir            VARCHAR(200) NOT NULL,
    UNIQUE (input_id),
    CONSTRAINT fk_pdbpredparams_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS af_msa;
CREATE TABLE af_msa (
    input_id   VARCHAR(20) COLLATE utf8mb4_bin NOT NULL,
    a3m_name   VARCHAR(60) NOT NULL,
    host       VARCHAR(10),
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    UNIQUE (input_id)
);

DROP TABLE IF EXISTS af_pred;
CREATE TABLE af_pred (
    input_id           VARCHAR(60) COLLATE utf8mb4_bin NOT NULL,
    ranking            INT,
    avg_plddt          FLOAT,
    avg_ptm            FLOAT,
    max_pae            FLOAT,
    json_name          VARCHAR(100),
    relaxed_pdb_name   VARCHAR(100),
    unrelaxed_pdb_name VARCHAR(100),
    UNIQUE (input_id, relaxed_pdb_name)
);





DROP TABLE IF EXISTS af_plddt;
CREATE TABLE af_plddt (
    input_id        VARCHAR(60) COLLATE utf8mb4_bin NOT NULL,
    relative_pdb    VARCHAR(300),
    relative_pkl    VARCHAR(300),
    model_no        INT,
    chain_no        VARCHAR(10) COLLATE utf8mb4_bin NOT NULL,
    avg_plddt       FLOAT,
    UNIQUE (input_id, relative_pdb, model_no, chain_no)
);

DROP TABLE IF EXISTS af_rmsd;
CREATE TABLE af_rmsd (
    input_id    VARCHAR(60) COLLATE utf8mb4_bin NOT NULL,
    pdb_status  VARCHAR(20),
    rmsd        FLOAT,
    native_pdb  VARCHAR(100),
    model_pdb   VARCHAR(300),
    UNIQUE(native_pdb, model_pdb)
);

DROP TABLE IF EXISTS af_tmalign;
CREATE TABLE af_tmalign(
    input_id    VARCHAR(60) COLLATE utf8mb4_bin NOT NULL,
    pdb_status  VARCHAR(20),
    rmsd        FLOAT,
    tm1         FLOAT,
    tm2         FLOAT,
    native_pdb  VARCHAR(100),
    model_pdb   VARCHAR(300),
    UNIQUE(native_pdb, model_pdb)
);
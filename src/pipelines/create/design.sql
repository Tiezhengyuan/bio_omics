/*
design_chain -> design_dimer  \
             -> design_trimer /
    
    -> predict -> predict_summary -> predict_plddt
                                -> unrelaxed_contacts, relaxed_contacts
                                -> unrelaxed_dist, relaxed_dist

*/

DROP TABLE IF EXISTS design_chain;
CREATE TABLE design_chain (
    target_name  VARCHAR(20),
    chain_id     VARCHAR(50) COLLATE utf8mb4_bin,
    parent_id    VARCHAR(30),
    pdb_id       VARCHAR(10),
    chain_no     VARCHAR(10) COLLATE utf8mb4_bin,
    faa_file     VARCHAR(100) NOT NULL,
    pdb_file     VARCHAR(200),
    binder_pdb   VARCHAR(200),
    chain_seq    TEXT,
    epitope_pos  VARCHAR(200),
    UNIQUE(target_name, chain_id)
);

DROP TABLE IF EXISTS design_dimer;
CREATE TABLE design_dimer (
    design_id       VARCHAR(100) COLLATE utf8mb4_bin PRIMARY KEY,
    category        VARCHAR(10),
    chain_labels    VARCHAR(20),
    target_name     VARCHAR(20),
    target_id       VARCHAR(40),
    target_chain_no VARCHAR(20) COLLATE utf8mb4_bin,
    binder_id       VARCHAR(40),
    binder_chain_no VARCHAR(20) COLLATE utf8mb4_bin,
    af_faa          VARCHAR(100) NOT NULL,
    target_pdb      VARCHAR(200),
    binder_pdb      VARCHAR(200),
    pdb_files       VARCHAR(1000)
);
DROP TABLE IF EXISTS design_trimer;
CREATE TABLE design_trimer LIKE design_dimer;

DROP TABLE IF EXISTS predict;
CREATE TABLE predict (
    predict_id     INT AUTO_INCREMENT PRIMARY KEY,
    target_name    VARCHAR(20),
    predictor      VARCHAR(10),
    design_id      VARCHAR(100),
    host           VARCHAR(10),
    endpoint       VARCHAR(100),
    start_at       TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    end_at         DATETIME,
    status         VARCHAR(10),
    UNIQUE(predictor, design_id)
);

DROP TABLE IF EXISTS predict_summary;
CREATE TABLE predict_summary (
    summary_id    VARCHAR(50) PRIMARY KEY,
    predict_id    INT REFERENCES predict(predict_id) ON DELETE CASCADE,
    ranking       INT,
    avg_plddt     FLOAT,
    avg_ptm       FLOAT,
    max_pae       FLOAT,
    relative_json VARCHAR(150),
    relaxed_pdb   VARCHAR(150),
    unrelaxed_pdb VARCHAR(150),
    UNIQUE (relaxed_pdb)
);

DROP TABLE IF EXISTS unrelaxed_plddt;
CREATE TABLE unrelaxed_plddt (
    summary_id    VARCHAR(50) REFERENCES predict_summary(summary_id) ON DELETE CASCADE,
    relative_pdb  VARCHAR(200),
    relative_pkl  VARCHAR(200),
    chain_no      VARCHAR(10),
    avg_plddt     FLOAT,
    UNIQUE (summary_id, chain_no)
);
DROP TABLE IF EXISTS relaxed_plddt;
CREATE TABLE relaxed_plddt LIKE unrelaxed_plddt;

DROP TABLE IF EXISTS unrelaxed_contacts;
CREATE TABLE unrelaxed_contacts(
    summary_id            VARCHAR(50) REFERENCES predict_summary(summary_id) ON DELETE CASCADE,
    relative_pdb          VARCHAR(200),
    binding_affinity      FLOAT,
    dissociation_constant FLOAT,
    relative_prodigy      VARCHAR(200),
    UNIQUE (summary_id)
);
DROP TABLE IF EXISTS relaxed_contacts;
CREATE TABLE relaxed_contacts LIKE unrelaxed_contacts;

DROP TABLE IF EXISTS unrelaxed_sasa;
CREATE TABLE unrelaxed_sasa(
    summary_id       VARCHAR(50) REFERENCES predict_summary(summary_id) ON DELETE CASCADE,
    delta_total_sasa FLOAT,
    relative_pkl     VARCHAR(200),
    relative_pdb     VARCHAR(200),
    UNIQUE (summary_id)
);
DROP TABLE IF EXISTS relaxed_sasa;
CREATE TABLE relaxed_sasa LIKE unrelaxed_sasa;

DROP TABLE IF EXISTS unrelaxed_dist;
CREATE TABLE unrelaxed_dist(
    summary_id    VARCHAR(50) REFERENCES predict_summary(summary_id) ON DELETE CASCADE,
    chain_combo   VARCHAR(10),
    min_dist      FLOAT,
    min_res1      VARCHAR(10),
    min_res2      VARCHAR(10),
    relative_pdb  VARCHAR(200),
    relative_pkl  VARCHAR(200),
    UNIQUE (summary_id, chain_combo)
);
DROP TABLE IF EXISTS relaxed_dist;
CREATE TABLE relaxed_dist LIKE unrelaxed_dist;

DROP TABLE IF EXISTS unrelaxed_dist_cdr;
CREATE TABLE unrelaxed_dist_cdr(
    summary_id  VARCHAR(50) REFERENCES predict_summary(summary_id) ON DELETE CASCADE,
    binder_id   VARCHAR(20),
    fragment    VARCHAR(10),
    res_no      INT,
    res         VARCHAR(10),
    min_dist    FLOAT
);
DROP TABLE IF EXISTS relaxed_dist_cdr;
CREATE TABLE relaxed_dist_cdr LIKE unrelaxed_dist_cdr;

-- DROP TABLE IF EXISTS docking;
-- CREATE TABLE docking (
--     docking_id     INT AUTO_INCREMENT PRIMARY KEY,
--     target_name    VARCHAR(20),
--     design_id      VARCHAR(100),
--     host           VARCHAR(10),
--     endpoint       VARCHAR(100),
--     start_at       TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
--     end_at         DATETIME,
--     status         VARCHAR(10),
--     UNIQUE(design_id)
-- );
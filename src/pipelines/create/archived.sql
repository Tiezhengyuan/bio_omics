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


/* drug protein */
DROP VIEW IF EXISTS view_drug;
CREATE VIEW view_drug AS
    SELECT C.chain_id, C.pdb_id,  C.uniprot_acc,
        C.chain_no, C.name, D.db_acc AS drug_acc,
        C.chain_seq, C.pro_seq, C.relative_faa, C.chain_pdb
    FROM view_pdb_uniprot     C
    LEFT JOIN uniprot U ON    C.uniprot_acc = U.uniprot_acc
    LEFT JOIN uniprot_dbxrefs D ON C.uniprot_acc = D.uniprot_acc
    WHERE U.organism LIKE '%Homo sapiens%'
        AND D.db_name = 'DrugBank'
;
//
DROP VIEW IF EXISTS view_ig;
CREATE VIEW view_ig AS
    SELECT I.allele_name, I.specie, I.gene_name, I.gene_family,
        I.chain_type, C.pdb_id, C.chain_id, C.first_chain_id,
        C.model_no, C.chain_no, C.pro_id, C.chain_seq, C.pro_len,
        C.release_date, C.structure_method, C.resolution, C.avg_bfactor,
        C.relative_faa, C.chain_pdb, C.pdb_group
    FROM view_pdb_chain     C
    LEFT JOIN align_vregion I ON C.chain_id = I.chain_id
    WHERE C.model_no = 0
        AND I.evalue <= 1e-5
        AND I.bit_score >= 50
        AND I.isotype = 'IG'
;
//
//
DROP VIEW IF EXISTS view_ab_dssp;
CREATE VIEW view_ab_dssp AS
    WITH tmp AS (
        SELECT I.allele_name, I.specie, I.gene_name, I.gene_family,
            I.chain_type, C.pdb_id, C.chain_id, C.model_no, C.chain_no,
            S.chain_seq, S.pro_len
        FROM align_vregion I, pdb_chainid C, pdb_chain_seq S
        WHERE I.chain_id = C.chain_id
            AND I.chain_id = S.chain_id
            AND I.evalue <= 1e-10
            AND I.bit_score >= 100
            AND I.isotype = 'IG'
    )
    SELECT A.pdb_id, A.chain_id, A.allele_name, A.specie,
        A.gene_name, A.gene_family, A.chain_type, B.avg_bfactor,
        A.chain_seq, B.bfactor, D.dssp_seq
    FROM tmp A
    LEFT JOIN pdb_bfactor B ON A.pdb_id=B.pdb_id 
        AND A.model_no=B.model_no AND A.chain_no=B.chain_no
    LEFT JOIN pdb_dssp D ON A.pdb_id=D.pdb_id
        AND A.model_no=D.model_no AND A.chain_no=D.chain_no
    WHERE B.bfactor is not null
        AND D.dssp_seq is not null
;
//
DROP VIEW IF EXISTS view_antibody_complex;
CREATE VIEW view_antibody_complex AS
    SELECT * FROM view_pdb_uniprot
    WHERE pdb_id IN (
        SELECT DISTINCT pdb_id FROM view_antibody
    )
    ORDER BY pdb_id
;
//
/* antibody sequences from ABSD aligned with IMGT V-region */
DROP VIEW IF EXISTS view_abseq_imgt;
CREATE VIEW view_abseq_imgt AS
    SELECT S.seq_id, S.specie, I.allele_name, I.gene_name, I.gene_family,
        I.chain_type, I.isotype, S.pro_seq, F.faa_file, A.aln_file
    FROM absd_seq S, absd_faa F, absd_igblastp A, absd_align_imgt I
    WHERE S.seq_id = F.seq_id
        AND S.seq_id = A.seq_id
        AND S.seq_id = I.seq_id
        AND I.evalue <= 1e-10
        AND I.bit_score >= 50
;
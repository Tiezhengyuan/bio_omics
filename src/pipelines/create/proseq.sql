/*
store protein sequences by split tables

proseq_100
proseq_110
proseq_120
proseq_400
proseq_1000
proseq_long

All column names could be identical
*/
DROP TABLE IF EXISTS proseq_100;
CREATE TABLE proseq_100 (
    seq_id       INT AUTO_INCREMENT PRIMARY KEY,
    source       VARCHAR(10),
    seq          VARCHAR(100) NOT NULL,
    seq_head     VARCHAR(10) GENERATED ALWAYS AS (SUBSTR(seq, 1, 10)) STORED,
    seq_len      INT GENERATED ALWAYS AS (LENGTH(seq)) STORED,
    relative_faa VARCHAR(60),
    UNIQUE(seq)
);

DROP TABLE IF EXISTS proseq_110;
CREATE TABLE proseq_110 (
    seq_id       INT AUTO_INCREMENT PRIMARY KEY,
    source       VARCHAR(10),
    seq          VARCHAR(110) NOT NULL,
    seq_head     VARCHAR(10) GENERATED ALWAYS AS (SUBSTR(seq, 1, 10)) STORED,
    seq_len      INT GENERATED ALWAYS AS (LENGTH(seq)) STORED,
    relative_faa VARCHAR(60),
    UNIQUE(seq)
);
CREATE INDEX index_proseq110 ON proseq_110 (seq_head);

DROP TABLE IF EXISTS proseq_120;
CREATE TABLE proseq_120 (
    seq_id       INT AUTO_INCREMENT PRIMARY KEY,
    source       VARCHAR(10),
    seq          VARCHAR(120) NOT NULL,
    seq_head     VARCHAR(10) GENERATED ALWAYS AS (SUBSTR(seq, 1, 10)) STORED,
    seq_len      INT GENERATED ALWAYS AS (LENGTH(seq)) STORED,
    relative_faa VARCHAR(60),
    UNIQUE(seq)
);
CREATE INDEX index_proseq120 ON proseq_120 (seq_head);

DROP TABLE IF EXISTS proseq_130;
CREATE TABLE proseq_130 (
    seq_id       INT AUTO_INCREMENT PRIMARY KEY,
    source       VARCHAR(10),
    seq          VARCHAR(130) NOT NULL,
    seq_head     VARCHAR(10) GENERATED ALWAYS AS (SUBSTR(seq, 1, 10)) STORED,
    seq_len      INT GENERATED ALWAYS AS (LENGTH(seq)) STORED,
    relative_faa VARCHAR(60),
    UNIQUE(seq)
);
CREATE INDEX index_proseq130 ON proseq_130 (seq_head);

DROP TABLE IF EXISTS proseq_400;
CREATE TABLE proseq_400 (
    seq_id       INT AUTO_INCREMENT PRIMARY KEY,
    source       VARCHAR(10),
    seq          VARCHAR(400) NOT NULL,
    seq_head     VARCHAR(10) GENERATED ALWAYS AS (SUBSTR(seq, 1, 10)) STORED,
    seq_len      INT GENERATED ALWAYS AS (LENGTH(seq)) STORED,
    relative_faa VARCHAR(60),
    UNIQUE(seq)
);
CREATE INDEX index_proseq400 ON proseq_400 (seq_head);

DROP TABLE IF EXISTS proseq_1000;
CREATE TABLE proseq_1000 (
    seq_id       INT AUTO_INCREMENT PRIMARY KEY,
    source       VARCHAR(10),
    seq          VARCHAR(1000) CHARACTER SET latin1,
    seq_head     VARCHAR(10) GENERATED ALWAYS AS (SUBSTR(seq, 1, 10)) STORED,
    seq_len      INT GENERATED ALWAYS AS (LENGTH(seq)) STORED,
    relative_faa VARCHAR(60),
    UNIQUE(seq)
);

DROP TABLE IF EXISTS proseq_long;
CREATE TABLE proseq_long (
    seq_id       INT AUTO_INCREMENT PRIMARY KEY,
    source       VARCHAR(10),
    seq          TEXT NOT NULL,
    seq_head     VARCHAR(10) GENERATED ALWAYS AS (SUBSTR(seq, 1, 10)) STORED,
    seq_len      INT GENERATED ALWAYS AS (LENGTH(seq)) STORED,
    relative_faa VARCHAR(60)
);
CREATE UNIQUE INDEX idx_proseqlong ON proseq_long (seq(255));

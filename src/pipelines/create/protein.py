#!/usr/bin/python
'''
requirements:
    - an empty database is created.
example: python src/pipelines/create_protein.py
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

SQL_TEXT = """
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
//
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
//
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
# //
# DROP TABLE IF EXISTS pro_pep;
# CREATE TABLE pro_pep (
#     pep_id      INT PRIMARY KEY,
#     seq         VARCHAR(10)
# );
# //
# DROP TABLE IF EXISTS protein_seq;
# CREATE TABLE protein_seq (
#     pro_id      INT PRIMARY KEY,
#     prefix      VARCHAR(3),
#     pro_len     INT,
#     pro_seq     TEXT
# );
# CREATE INDEX index_prolen on protein_seq (pro_len);
# //
# DROP TABLE IF EXISTS protein_source;
# CREATE TABLE protein_source (
#     pro_id          INT,
#     source_name     VARCHAR(50),
#     source_id       VARCHAR(50),
#     UNIQUE(source_name, source_id),
#     CONSTRAINT fk_protein_source FOREIGN KEY (pro_id)
#         REFERENCES protein_seq(pro_id) ON DELETE CASCADE
# );
"""

#####################################################################################

if __name__ == "__main__":
    bc = BuildComplex(verbose=True)
    for query in SQL_TEXT.split('//'): 
        res = bc.create_table(query)
        if res is None:
            print(f"ERROR. Check the above SQL: {query}\n\n")

#!/usr/bin/python
'''
example: python src/pipelines/db_complex_create.py
functions:
    - build database known as complex2
    - build tables
'''
from Bio import SeqIO
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)


from bioomics import BuildComplex

# PDB data
create_pdb = """
    CREATE TABLE pdb(
        pdb_id VARCHAR(10),
        pdb_url VARCHAR(50),
        local_path VARCHAR(100),
        release_date DATE,
        resolution float,
        name VARCHAR(500),
        CONSTRAINT pk_pdb PRIMARY KEY (pdb_id)
    );
"""

create_compound = """
    CREATE TABLE compound(
        pdb_id VARCHAR(10),
        compound_no INT,
        chain VARCHAR(200),
        compound_text Text,
        CONSTRAINT pk_compound PRIMARY KEY (pdb_id, compound_no),
        CONSTRAINT fk_compound_pdb FOREIGN KEY (pdb_id)
            REFERENCES pdb(pdb_id)  ON DELETE CASCADE
    );
"""

create_chain = """
    CREATE TABLE chain(
        pdb_id VARCHAR(10),
        compound_no INT,
        chain_id VARCHAR(10),
        pro_seq Text,
        CONSTRAINT pk_chain PRIMARY KEY (pdb_id, compound_no, chain_id),
        CONSTRAINT fk_chain_pdb FOREIGN KEY (pdb_id)
            REFERENCES pdb(pdb_id)  ON DELETE CASCADE,
        CONSTRAINT fk_chain_compound FOREIGN KEY (pdb_id, compound_no)
            REFERENCES compound(pdb_id, compound_no)  ON DELETE CASCADE
    );
"""

create_pdb_keyword = """
    CREATE TABLE pdb_keyword(
        pdb_id VARCHAR(10),
        keyword VARCHAR(50),
        UNIQUE (pdb_id, keyword),
        CONSTRAINT fk_pdb_keyword_pdb FOREIGN KEY (pdb_id)
            REFERENCES pdb(pdb_id) ON DELETE CASCADE
    );
"""

# uniprot data: accession primary key
create_uniprot = """
    CREATE TABLE uniprot (
        uniprot_acc VARCHAR(20),
        name VARCHAR(100),
        organism VARCHAR(200),
        pro_seq TEXT,
        CONSTRAINT pk_uniprot PRIMARY KEY (uniprot_acc)
    );
"""

# uniprot data: keywords for aggregating data
create_uniprot_keyword = """
    CREATE TABLE uniprot_keyword (
        uniprot_acc VARCHAR(20),
        keyword VARCHAR(100),
        UNIQUE (uniprot_acc, keyword),
        CONSTRAINT fk_uniprot_keyword_uniprot FOREIGN KEY (uniprot_acc)
            REFERENCES uniprot(uniprot_acc) ON DELETE CASCADE
    );
"""

# uniprot data: set links between uniprot and other database
create_uniprot_dbxrefs = """
    CREATE TABLE uniprot_dbxrefs (
        uniprot_acc VARCHAR(20),
        db_name VARCHAR(50),
        db_acc VARCHAR(50),
        UNIQUE (uniprot_acc, db_name, db_acc),
        CONSTRAINT fk_uniprot_dbxrefs_uniprot FOREIGN KEY (uniprot_acc)
            REFERENCES uniprot(uniprot_acc) ON DELETE CASCADE
    );
"""

# blastp of pdb chains on uniprot_sprot
create_uniprot_align = """
    CREATE TABLE uniprot_align(
        id          INT AUTO_INCREMENT PRIMARY KEY,
        pdb_id      VARCHAR(10),
        model_id    INT,
        chain_id    VARCHAR(5),
        uniprot_acc VARCHAR(20),
        evalue      DOUBLE,
        bit_score   INT,
        identity    FLOAT
    );
"""
# views
view_pdb_uniprot = """
    CREATE VIEW view_pdb_uniprot AS
        SELECT A.pdb_id, A.chain_id, B.uniprot_acc, B.name, C.keyword, D.pro_seq
        FROM uniprot_align A
        LEFT JOIN uniprot B
        ON A.uniprot_acc = B.uniprot_acc
        LEFT JOIN uniprot_keyword C
        ON A.uniprot_acc = C.uniprot_acc
        LEFT JOIN chain D
        on A.pdb_id = D.pdb_id AND A.chain_id = D.chain_id
        WHERE A.model_id = 0 AND bit_score > 100 AND evalue < 1e-10
            AND D.pro_seq IS NOT NULL;
"""

view_antibody = """
    CREATE VIEW view_antibody AS
        SELECT A.pdb_id, A.chain_id, B.uniprot_acc, B.name, C.keyword, D.pro_seq
        FROM uniprot_align A
        LEFT JOIN uniprot B
        ON A.uniprot_acc = B.uniprot_acc
        LEFT JOIN uniprot_keyword C
        ON A.uniprot_acc = C.uniprot_acc
        LEFT JOIN chain D
        on A.pdb_id = D.pdb_id AND A.chain_id = D.chain_id
        WHERE A.model_id = 0 AND A.bit_score > 100 AND A.evalue < 1e-10
            AND C.keyword REGEXP 'immunoglobulin|antibody'
            AND D.pro_seq IS NOT NULL;
"""

view_antibody_complex = """
    CREATE VIEW view_antibody_complex AS
        SELECT DISTINCT pdb_id, chain_id, name, pro_seq
        FROM view_pdb_uniprot
        where pdb_id IN (
            SELECT DISTINCT pdb_id
            FROM view_antibody
        );
"""

# igblastp of pdb chains on V-gene
create_abag = """
    CREATE TABLE abag (
        pdb_id VARCHAR(10),
        compound_no INT,
        chain_id VARCHAR(10),
        chain_type VARCHAR(30),
        fragment VARCHAR(30),
        CONSTRAINT pk_abag PRIMARY KEY (pdb_id, chain_id),
        CONSTRAINT fk_abag_compound FOREIGN KEY (pdb_id, compound_no)
            REFERENCES compound(pdb_id, compound_no) ON DELETE CASCADE
    );
"""

create_imgt_genelist = """
    CREATE TABLE imgt_genelist (
        species VARCHAR(50),
        gene_db VARCHAR(20),
        gene_function VARCHAR(20),
        gene_definition VARCHAR(200),
        number_alleles VARCHAR(5),
        chromosome VARCHAR(5),
        chromosomal_loc VARCHAR(20),
        ligm_db_ref_gene VARCHAR(100),
        ncbi_gene VARCHAR(20),
        chain VARCHAR(10),
        domain VARCHAR(5)
    );
"""



#####################################################################################
def initialize_tables(meta):
    b = BuildComplex(verbose=True)

    # drop tables
    for table_name, create_query in meta[::-1]:
        res = b.drop_table(table_name)
        print(f'Try to drop table {table_name}:', res)
    b.drop_database()

    # create tables
    b.create_database()
    for table_name, create_query in meta:
        res = b.create_table(create_query)
        print(f'Try to create table {table_name}:', res)

if __name__ == "__main__":
    meta = [
        # pdb determined
        ('pdb', create_pdb),
        ('pdb_keyword', create_pdb_keyword),
        ('compound', create_compound),
        ('chain', create_chain),
        # pdb related
        ('abag', create_abag),
        ('uniprot_align', create_uniprot_align),
        # other db
        ('uniprot', create_uniprot),
        ('uniprot_keyword', create_uniprot_keyword),
        ('uniprot_dbxrefs', create_uniprot_dbxrefs),
        # view
        ('view_pdb_uniprot', view_pdb_uniprot),
    ]
    initialize_tables(meta)

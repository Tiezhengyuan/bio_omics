/*
view_proseq
view_pdb -> view_pdb_chain -> view_pdb_uniprot
                           -> view_complex
                            -> view_antibody
*/
DROP VIEW IF EXISTS view_proseq;
CREATE VIEW view_proseq AS
    SELECT 's100' AS proseq_name, seq_id, source, seq_head, seq_len, seq, relative_faa
    FROM proseq_400
    UNION
    SELECT 's110' AS proseq_name, seq_id, source, seq_head, seq_len, seq, relative_faa
    FROM proseq_400
    UNION
    SELECT 's120' AS proseq_name, seq_id, source, seq_head, seq_len, seq, relative_faa
    FROM proseq_400
    UNION
    SELECT 's130' AS proseq_name, seq_id, source, seq_head, seq_len, seq, relative_faa
    FROM proseq_400
    UNION
    SELECT 's400' AS proseq_name, seq_id, source, seq_head, seq_len, seq, relative_faa
    FROM proseq_400
    UNION
    SELECT 's1000' AS proseq_name, seq_id, source, seq_head, seq_len, seq, relative_faa
    FROM proseq_1000
    UNION
    SELECT 'long' AS proseq_name, seq_id, source, seq_head, seq_len, seq, relative_faa
    FROM proseq_long
;

DROP VIEW IF EXISTS view_pdb;
CREATE VIEW view_pdb AS
    SELECT P.pdb_id, P.pdb_url, B.relative_pdb, B.pdb_group,
        P.release_date, P.name, M.structure_method, M.resolution
    FROM pdb P, meta_pdb B, pdb_method M
    WHERE P.pdb_id = B.pdb_id
        AND P.pdb_id = M.pdb_id
        AND B.relative_pdb IS NOT NULL
        AND M.resolution IS NOT NULL
;

DROP VIEW IF EXISTS view_pdb_chain;
CREATE VIEW view_pdb_chain AS
    SELECT C.chain_id, C.pdb_id, C.model_no, C.chain_no,
        S.first_chain_id, X.proseq_name, X.seq_id,
        D.structure_method, D.resolution, B.avg_bfactor,
        D.release_date, D.pdb_group,
        S.chain_seq, S.pro_len,
        S.relative_faa, P.relative_pdb
    FROM pdb_chainid        C
    LEFT JOIN pdb_chain_seq S ON C.chain_id = S.chain_id
    LEFT JOIN chain_pdb     P ON C.chain_id = P.chain_id
    LEFT JOIN view_pdb      D ON C.pdb_id = D.pdb_id
    LEFT JOIN pdb_proseq    X ON C.chain_id = X.chain_id
    LEFT JOIN pdb_bfactor   B ON C.pdb_id = B.pdb_id 
        AND C.model_no = B.model_no AND C.chain_no = B.chain_no
    WHERE S.chain_seq IS NOT NULL
;

DROP VIEW IF EXISTS view_complex;
CREATE VIEW view_complex AS
    SELECT pdb_id, count(*) AS num_chains, GROUP_CONCAT(DISTINCT 
        chain_type ORDER BY chain_type) AS complex_type
    FROM pdb_chain_type
    GROUP BY pdb_id
    ORDER BY num_chains, complex_type
;

/* antibody */
DROP VIEW IF EXISTS view_antibody;
CREATE VIEW view_antibody AS
    SELECT I.allele_name, I.specie, I.gene_name, I.gene_family,
        I.chain_type, C.pdb_id, C.chain_id, C.first_chain_id,
        C.model_no, C.chain_no, C.proseq_name, C.seq_id, C.chain_seq,
        C.pro_len, C.release_date, C.structure_method, C.resolution,
        C.avg_bfactor, C.relative_faa, C.relative_pdb, C.pdb_group
    FROM view_pdb_chain     C
    LEFT JOIN align_vregion I ON C.chain_id = I.chain_id
    WHERE C.model_no = 0
        AND I.evalue <= 1e-10
        AND I.bit_score >= 100
        AND I.isotype = 'IG'
        AND I.region_names = 'FR1,CDR1,FR2,CDR2,FR3,CDR3'
;

-- T-cell receptor
-- DROP VIEW IF EXISTS view_tcr;
-- CREATE VIEW view_tcr AS
--     SELECT I.allele_name, I.gene_name, I.gene_family,
--         I.chain_type, C.pdb_id, C.chain_id, C.model_no,
--         C.chain_no, C.seq_id, C.chain_seq, C.structure_method,
--         C.resolution, C.relative_faa, C.relative_pdb, C.pdb_group
--     FROM view_pdb_chain C, align_vregion I
--     WHERE C.chain_id = I.chain_id
--         AND I.evalue <= 1e-10
--         AND I.bit_score >= 50
--         AND I.isotype = 'TR'
-- ;   

/* uniprot */
DROP VIEW IF EXISTS view_pdb_uniprot;
CREATE VIEW view_pdb_uniprot AS
    SELECT A.chain_id, A.pdb_id, A.uniprot_acc, U.pro_name,
        C.chain_no, C.chain_seq,
        A.identity, U.pro_seq,
        C.structure_method, C.resolution,
        C.relative_faa, C.relative_pdb
    FROM align_uniprot       A
    LEFT JOIN uniprot        U ON A.uniprot_acc = U.uniprot_acc
    LEFT JOIN view_pdb_chain C ON A.chain_id = C.chain_id
    WHERE A.bit_score >= 50
        AND A.evalue <= 1e-10
        AND A.identity >= .3
        AND C.model_no = 0
;

/* human transmembrane proteins */
-- DROP VIEW IF EXISTS view_htrm;
-- CREATE VIEW view_htrm AS
--     SELECT C.chain_id, C.pdb_id,  C.uniprot_acc,
--         C.chain_no, C.pro_name, C.chain_seq, C.pro_seq
--     FROM view_pdb_uniprot C
--     LEFT JOIN uniprot U
--     ON C.uniprot_acc = U.uniprot_acc
--     WHERE U.organism LIKE '%Homo sapiens%'
--         AND C.uniprot_acc IN (
--             SELECT distinct uniprot_acc
--             FROM uniprot_keyword
--             WHERE keyword LIKE 'Transmembrane'
--         )
-- ;

/* allergen */
-- DROP VIEW IF EXISTS view_allergen;
-- CREATE VIEW view_allergen AS
--     SELECT C.chain_id, C.pdb_id,  C.uniprot_acc,
--         C.chain_no, C.pro_name, D.db_acc AS allergen_acc,
--         C.chain_seq, C.pro_seq, C.relative_faa , C.relative_pdb
--     FROM view_pdb_uniprot     C
--     LEFT JOIN uniprot         U ON C.uniprot_acc = U.uniprot_acc
--     LEFT JOIN uniprot_dbxrefs D ON C.uniprot_acc = D.uniprot_acc
--     WHERE D.db_name = 'allergome'
--         AND U.organism LIKE '%Homo sapiens%'
-- ;
#! /usr/bin/bash

APP="uv --directory /home/yuan/bio/bio_omics run bomics"


######### alphafold
# build params
parallel -j8 ${APP} af_params ::: {0..241}
# MSA
parallel -j8 ${APP} af_msa ::: {0..241}
# predict
parallel -j8 ${APP} af_pred ::: {0..241}
# summary
parallel -j8 ${APP} af_summary ::: {0..241}
# pLDDT
parallel -j8 ${APP} af_plddt ::: {0..241}
# RMSD
parallel -j8 ${APP} af_rmsd ::: {0..241}
# TM-align score
parallel -j8 ${APP} af_tmalign ::: {0..241}



SELECT a.combo_id AS input_id,
    STRING_AGG(a.category, ',')
FROM label_complex a
LEFT JOIN meta_pdb      p ON a.pdb_id = p.pdb_id
WHERE p.pdb_group = 0
    AND a.combo = 5
    AND NOT EXISTS (
        SELECT 1 FROM af_params b
        WHERE a.combo_id = b.input_id
    )
GROUP BY a.combo_id;



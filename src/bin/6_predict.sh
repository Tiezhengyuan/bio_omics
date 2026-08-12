#! /usr/bin/bash

APP="uv --directory /home/yuan/bio/bio_omics run bomics"

parallel -j8 ${APP} pdb_pred_params ::: {0..241}


######### alphafold
# MSA
parallel -j8 ${APP} af_msa ::: {0..241}
# predict
parallel -j8 ${APP} af_pred ::: {0..241}
# pLDDT
parallel -j8 ${APP} af_plddt ::: {0..241}
# RMSD
parallel -j8 ${APP} af_rmsd ::: {0..241}
# TM-align score
parallel -j8 ${APP} af_tmalign ::: {0..241}


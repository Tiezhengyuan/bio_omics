#! /usr/bin/bash

APP=/home/yuan/bio/bio_omics/app.py

# alphafold
parallel -j8 python ${APP} chain_alphafold ::: {0..117}
parallel -j8 python ${APP} chain_alphafold_plddt ::: {0..117}
parallel -j8 python ${APP} chain_afsum ::: {0..117}
parallel -j8 python ${APP} chain_eval ::: {0..117}

parallel -j8 python ${APP} combo2_alphafold ::: {0..117}
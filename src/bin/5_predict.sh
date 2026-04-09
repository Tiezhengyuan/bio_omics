#! /usr/bin/bash

APP=/home/yuan/bio/bio_omics/app.py

###monomer
# alphafold
parallel -j8 python ${APP} chain_af ::: {0..237}
parallel -j8 python ${APP} chain_af_summary ::: {0..237}
# update table chain_plddt
parallel -j8 python ${APP} chain_af_plddt ::: {0..237}

# need chain_plddt
# update table chain_rmsd
parallel -j8 python ${APP} chain_rmsd ::: {0..237}---
# updat table chain tmalign
parallel -j8 python ${APP} chain_tmalign ::: {0..237}---


###dimer
parallel -j8 python ${APP} combo2_af ::: {0..237}---

parallel -j8 python ${APP} combo3_af ::: {0..237}---

parallel -j8 python ${APP} combo4_af ::: {0..237}---

parallel -j8 python ${APP} combo5_af ::: {0..237}---
#! /usr/bin/bash

APP=/home/yuan/bio/bio_omics/app.py


###monomer: need label_chain, chain_relaxed
# predict chain using alphafold
parallel -j8 python ${APP} chain_af ::: {0..240}
parallel -j8 python ${APP} chain_af_summary ::: {0..240}
# update chain_plddt
parallel -j8 python ${APP} chain_af_plddt ::: {0..240}
# need chain_af_plddt
parallel -j8 python ${APP} chain_rmsd ::: {0..240}
parallel -j8 python ${APP} chain_tmalign ::: {0..240}

### multimer: need label_complex
#dimer
parallel -j8 python ${APP} combo2_af ::: {0..240}---

#trimer
parallel -j8 python ${APP} combo3_af ::: {0..240}---

#tetramer
parallel -j8 python ${APP} combo4_af ::: {0..240}---

#pentamer
parallel -j8 python ${APP} combo5_af ::: {0..240}---
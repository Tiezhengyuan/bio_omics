#! /usr/bin/bash

APP="uv --directory /home/yuan/bio/bio_omics run bomics"


###monomer: need label_chain, chain_relaxed
# predict chain using alphafold
parallel -j8 ${APP} chain_af ::: {0..241}
parallel -j8 ${APP} chain_af_summary ::: {0..241}
# update chain_plddt
parallel -j8 ${APP} chain_af_plddt ::: {0..241}
# need chain_af_plddt
parallel -j8 ${APP} chain_rmsd ::: {0..241}
parallel -j8 ${APP} chain_tmalign ::: {0..241}


### multimer: need label_complex
#dimer
parallel -j8 ${APP} combo2_af ::: {0..241}---

#trimer
parallel -j8 ${APP} combo3_af ::: {0..241}---

#tetramer
parallel -j8 ${APP} combo4_af ::: {0..241}---

#pentamer
parallel -j8 ${APP} combo5_af ::: {0..241}---
#! /usr/bin/bash

APP="uv --directory /home/yuan/bio/bio_omics run bomics"


# need label_chain, meta_chain_pdb
parallel -j8 ${APP} chain_fixed ::: {0..241}
parallel -j8 ${APP} chain_relaxed ::: {0..241}

# RMSD of various structure with identical seq
parallel -j8 ${APP} app.py equal_seq_rmsd ::: {0..241}

#######NOTE: label_complex should be ready

# dimer
parallel -j8 ${APP} combo2 ::: {0..241}
parallel -j8 ${APP} combo2_fixed ::: {0..241}
parallel -j8 ${APP} combo2_relaxed ::: {0..241}

# trimer
parallel -j8 ${APP} combo3 ::: {0..241}
parallel -j8 ${APP} combo3_fixed ::: {0..241}
parallel -j8 ${APP} combo3_relaxed ::: {0..241}

#tetramer
parallel -j8 ${APP} combo4 ::: {0..241}
parallel -j8 ${APP} combo4_fixed ::: {0..241}
parallel -j8 ${APP} combo4_relaxed ::: {0..241}

#pentamer
parallel -j8 ${APP} combo5 ::: {0..241}
parallel -j8 ${APP} combo5_fixed ::: {0..241}
parallel -j8 ${APP} combo5_relaxed ::: {0..241}

#geometry 
parallel -j8 ${APP} pdb_qc ::: {0..241}
parallel -j8 ${APP} pdb_relaxed_qc ::: {0..241}
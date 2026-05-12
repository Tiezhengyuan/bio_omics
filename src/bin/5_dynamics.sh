#! /usr/bin/bash

APP="/home/yuan/bio/bio_omics/app.py"


# need label_chain, meta_chain_pdb
parallel -j8 python ${APP} chain_relaxed ::: {0..240}

# RMSD of various structure with identical seq
parallel -j8 python ${APP} app.py equal_seq_rmsd ::: {0..240}

#######NOTE: label_complex should be ready

# dimer
parallel -j8 python ${APP} combo2 ::: {0..240}
# need meta_combo2
parallel -j8 python ${APP} combo2_relaxed ::: {0..240}

# trimer
parallel -j8 python ${APP} combo3 ::: {0..240}
# need meta_combo3
parallel -j8 python ${APP} combo3_relaxed ::: {0..240}

#tetramer
parallel -j8 python ${APP} combo4 ::: {0..240}
# need meta_combo4
parallel -j8 python ${APP} combo4_relaxed ::: {0..240}

#pentamer
parallel -j8 python ${APP} combo5 ::: {0..240}
# need meta_combo5
parallel -j8 python ${APP} combo5_relaxed ::: {0..240}


# atom distance
parallel -j8 python ${APP} combo2_distance ::: {0..240}

# SASA
parallel -j8 python ${APP} combo_sasa ::: {0..240}

# simulate binding affinity
parallel -j8 python ${APP} combo_contacts ::: {0..240}

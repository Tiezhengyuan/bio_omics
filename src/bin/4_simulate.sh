#! /usr/bin/bash

APP="/home/yuan/bio/bio_omics/app.py"

#######NOTE: meta_pdb, view_antibody, pdb_complex should be ready

# chain structures
parallel -j8 python ${APP} meta_chain_pdb ::: {0..237}
# need meta_chain_pdb
parallel -j8 python ${APP} meta_chain_relaxed ::: {0..237}

# dimer
parallel -j8 python ${APP} combo2 ::: {0..237}
# need meta_combo2
parallel -j8 python ${APP} combo2_relaxed ::: {0..237}---

# trimer
parallel -j8 python ${APP} combo3 ::: {0..237}
# need meta_combo3
parallel -j8 python ${APP} combo3_relaxed ::: {0..237}---

#tetramer
parallel -j8 python ${APP} combo4 ::: {0..237}
# need meta_combo4
parallel -j8 python ${APP} combo4_relaxed ::: {0..237}

#pentamer
parallel -j8 python ${APP} combo5 ::: {0..237}
# need meta_combo5
parallel -j8 python ${APP} combo5_relaxed ::: {0..237}



# meta_combo2
# atom distance
parallel -j8 python ${APP} combo2_distance ::: {0..237}---
# SASA
parallel -j8 python ${APP} combo2_sasa ::: {0..237}---



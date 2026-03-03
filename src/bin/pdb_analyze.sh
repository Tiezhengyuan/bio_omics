#! /usr/bin/bash

APP_DIR=/home/yuan/bio/bio_omics

# calculate B-factor
parallel -j8 python ${APP_DIR}/app.py meta_bfactor ::: {0..117}
parallel -j8 python ${APP_DIR}/app.py pdb_bfactor ::: {0..117}

# identify chain type
parallel -j8 python ${APP_DIR}/app.py pdb_chain_type ::: {0..117}

# secondary structure
parallel -j8 python ${APP_DIR}/app.py pdb_dssp ::: {0..117}

# atom distance
parallel -j8 python ${APP_DIR}/app.py combo2_distance ::: {0..117}
# require meta_combo_pdb
parallel -j8 python ${APP_DIR}/app.py pdb_sasa ::: {0..117}

# RMSD
parallel -j8 python ${APP_DIR}/app.py equal_seq_rmsd ::: {0..117}

#contacts
parallel -j8 python ${APP_DIR}/app.py combo2_contacts ::: {0..117}

# antibody numbering
parallel -j8 python ${APP_DIR}/app.py anarci ::: {0..117}
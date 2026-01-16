#! /usr/bin/bash

PIPELINE_DIR=/home/yuan/bio/bio_omics/src/pipelines

# calculate B-factor
parallel -j8 python ${PIPELINE_DIR}/meta_bfactor.py ::: {0..117}
parallel -j8 python ${PIPELINE_DIR}/pdb_bfactor.py ::: {0..117}

# identify chain type
parallel -j8 python ${PIPELINE_DIR}/pdb_chain_type.py ::: {0..117}

# secondary structure
parallel -j8 python ${PIPELINE_DIR}/pdb_dssp.py ::: {0..117}

# atom distance
parallel -j8 python ${PIPELINE_DIR}/combo2_distance.py ::: {0..117}
# require meta_combo_pdb
parallel -j8 python ${PIPELINE_DIR}/pdb_sasa.py ::: {0..117}

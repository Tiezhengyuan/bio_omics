#! /usr/bin/bash

PIPELINE_DIR=/home/yuan/bio/bio_omics/src/pipelines

# calculate B-factor
parallel -j4 python ${PIPELINE_DIR}/meta_bfactor.py ::: {0..116}
parallel -j8 python ${PIPELINE_DIR}/pdb_bfactor.py ::: {0..116}

# identify chain type
parallel -j8 python ${PIPELINE_DIR}/pdb_chain_type.py ::: {0..116}

# secondary structure
parallel -j4 python ${PIPELINE_DIR}/pdb_dssp.py ::: {0..116}

# atom distance
parallel -j4 python ${PIPELINE_DIR}/pdb_distance.py ::: {0..116}
# require meta_combo_pdb
parallel -j4 python ${PIPELINE_DIR}/pdb_sasa.py ::: {0..116}



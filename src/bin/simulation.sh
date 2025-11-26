#! /usr/bin/bash

PIPELINE_DIR=/home/yuan/bio/bio_omics/src/pipelines

# fix pdb
parallel -j4 python ${PIPELINE_DIR}/pdb_fixer.py ::: {0..46}

parallel -j4 python ${PIPELINE_DIR}/pdb_fixed.py ::: {0..46}
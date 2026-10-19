#! /usr/bin/bash

APP="uv --directory /home/yuan/bio/bio_omics run bomics"

# RMSD of various structure with identical seq
parallel -j8 ${APP} equal_seq_rmsd ::: {0..241}

# NOTE: chain_fixed, chain_relaxed, combo2..5 and their *_fixed/*_relaxed
# pipelines moved to structure_profiling:
#   see ~/bio/structure_profiling/bin/5_dynamics.sh

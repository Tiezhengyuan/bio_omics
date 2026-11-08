#! /usr/bin/bash

APP="uv --directory /home/yuan/bio/bio_omics run bomics"

# NOTE: equal_seq_rmsd, chain_fixed, chain_relaxed, combo2..5 and their
# *_fixed/*_relaxed pipelines moved to structure_profiling:
#   equal_seq_rmsd:
#     ~/bio/structure_profiling/src/pipelines/dynamics/equal_seq_rmsd.py
#     run: sp equal_seq_rmsd (see ~/bio/structure_profiling/bin/5_qc.sh)

#! /usr/bin/bash
# example: bash src/bin/pdb_relaxe.sh /home/yuan/bio/bio_omics/src/pipelines

PIPELINE_DIR=$1

# fix pdb
parallel -j4 python ${PIPELINE_DIR}/fixed_pdb.py ::: {0..116}
# split fixed_pdb
parallel -j4 python ${PIPELINE_DIR}/fixed_chainid.py ::: {0..116}
parallel -j4 python ${PIPELINE_DIR}/fixed_chain.py ::: {0..116}
parallel -j4 python ${PIPELINE_DIR}/fixed_combo2.py ::: {0..116}
parallel -j4 python ${PIPELINE_DIR}/fixed_combo3.py ::: {0..116}

# from chain 
parallel -j4 python ${PIPELINE_DIR}/pdb_chain_fixed.py ::: {0..116}

# relax pdb then split PDB
parallel -j4 python ${PIPELINE_DIR}/relaxed_pdb_amber99sb.py ::: {0..116}
parallel -j4 python ${PIPELINE_DIR}/relaxed_pdb_amber14.py ::: {0..116}
parallel -j4 python ${PIPELINE_DIR}/relaxed_amber14_chain.py ::: {0..116}

# relax chain
parallel -j4 python ${PIPELINE_DIR}/chain_relaxed_amber14.py ::: {0..116}

# fix + relax
parallel -j4 python ${PIPELINE_DIR}/meta_chain_relaxed.py ::: {0..116}
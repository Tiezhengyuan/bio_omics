#! /usr/bin/bash

PIPELINE_DIR=/home/yuan/bio/bio_omics/src/pipelines

# synchronize raw data from PDB
python ${PIPELINE_DIR}/pdb_rsync.py

# PDB to pdb
python ${PIPELINE_DIR}/put_pdb.py
# pdb_files of PDB
python ${PIPELINE_DIR}/put_meta_pdb.py

# require meta_pdb
parallel -j8 python ${PIPELINE_DIR}/put_compound.py ::: {0..46}
parallel -j8 python ${PIPELINE_DIR}/put_method.py ::: {0..46}
parallel -j8 python ${PIPELINE_DIR}/put_keyword.py ::: {0..46}
parallel -j8 python ${PIPELINE_DIR}/put_chain_type.py ::: {0..46}
parallel -j8 python ${PIPELINE_DIR}/put_chain.py ::: {0..46}
parallel -j8 python ${PIPELINE_DIR}/put_chainid.py ::: {0..46}

# require pdb_chain
parallel -j8 python ${PIPELINE_DIR}/put_chain_seq.py ::: {0..46}
parallel -j8 python ${PIPELINE_DIR}/export_faa.py ::: {0..46}
parallel -j8 python ${PIPELINE_DIR}/export_chain_faa.py ::: {0..46}
parallel -j8 python ${PIPELINE_DIR}/export_chain_pdb.py ::: {0..46}
parallel -j8 python ${PIPELINE_DIR}/export_pair_pdb.py ::: {0..46}

# pdb_chain should be ready
# calculate B-factor
parallel -j4 python ${PIPELINE_DIR}/pdb_bfactor.py ::: {0..46}
parallel -j8 python ${PIPELINE_DIR}/put_bfactor.py ::: {0..46}

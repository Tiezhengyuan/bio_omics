#! /usr/bin/bash

PIPELINE_DIR=/home/yuan/bio/bio_omics/src/pipelines

# synchronize raw data from PDB
python ${PIPELINE_DIR}/download_pdb.py
# python ${PIPELINE_DIR}/download_mmcif.py

# PDB to pdb
python ${PIPELINE_DIR}/put_pdb.py 
# python ${PIPELINE_DIR}/put_pdb2.py

# pdb_files of PDB
python ${PIPELINE_DIR}/put_meta_pdb.py
# python ${PIPELINE_DIR}/put_meta_mmcif.py

# require meta_pdb
parallel -j8 python ${PIPELINE_DIR}/pdb_compound.py ::: {0..116}
parallel -j8 python ${PIPELINE_DIR}/pdb_method.py ::: {0..116}
parallel -j8 python ${PIPELINE_DIR}/pdb_keyword.py ::: {0..116}
# parallel -j8 python ${PIPELINE_DIR}/put_chain.py ::: {0..116}
parallel -j8 python ${PIPELINE_DIR}/pdb_chainid.py ::: {0..116}

# require pdb_chain
parallel -j8 python ${PIPELINE_DIR}/pdb_chain_seq.py ::: {0..116}
parallel -j8 python ${PIPELINE_DIR}/unique_chain_seq.py
parallel -j8 python ${PIPELINE_DIR}/meta_faa.py ::: {0..116}
parallel -j8 python ${PIPELINE_DIR}/meta_chain_faa.py ::: {0..116}
parallel -j8 python ${PIPELINE_DIR}/meta_chain_pdb.py ::: {0..116}

parallel -j8 python ${PIPELINE_DIR}/meta_combo2_pdb.py ::: {0..116}
# file for prediciton: .csv for alphafold
parallel -j8 python ${PIPELINE_DIR}/meta_combo2_seq.py ::: {0..116}

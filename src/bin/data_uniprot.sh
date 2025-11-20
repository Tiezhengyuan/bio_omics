#! /usr/bin/bash

PIPELINE_DIR=/home/yuan/bio/bio_omics/src/pipelines

# load uniprot_acc
python ${PIPELINE_DIR}/uniprot.py

# load uniport annotations
python ${PIPELINE_DIR}/uniprot_keyword.py &
P1=$!
python ${PIPELINE_DIR}/uniprot_dbxrefs.py &
P2=$!
python ${PIPELINE_DIR}/uniprot_feature.py &
P3=$!
python ${PIPELINE_DIR}/uniprot_meta_alphafold.py &
P4=$!
wait $P1 $P2 $P3 $P4

# pLDDT of alphafold modelling
python ${PIPELINE_DIR}/uniprot_alphafold_plddt.py

# RMSD between chain pdb and uniprot prediction by alphafold
parallel -j4 python ${PIPELINE_DIR}/uniprot_alphafold_rmsd.py ::: {0..46}
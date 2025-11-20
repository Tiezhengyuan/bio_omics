#! /usr/bin/bash
echo "Try to integrate unique protein sequences..."
PIPELINE_DIR=/home/yuan/bio/bio_omics/src/pipelines

echo 'Data from PDB...'
python ${PIPELINE_DIR}/protein_pdb.py
echo 'Data from Uniprot...'
python ${PIPELINE_DIR}/protein_uniprot.py
echo 'Data from IMGT/INN...'
python ${PIPELINE_DIR}/protein_inn.py


#! /usr/bin/bash

PIPELINE_DIR=/home/yuan/bio/bio_omics/src/pipelines
LIST_DIR=/home/yuan/data/pdb/_list


python ${PIPELINE_DIR}/uniprot.py

python ${PIPELINE_DIR}/uniprot_keyword.py &
P1=$!
python ${PIPELINE_DIR}/uniprot_dbxrefs.py &
P2=$!
python ${PIPELINE_DIR}/uniprot_meta_alphafold.py &
P3=$!
wait $P1 $P2 $P3

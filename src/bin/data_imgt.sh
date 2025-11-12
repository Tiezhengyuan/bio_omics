#! /usr/bin/bash

PIPELINE_DIR=/home/yuan/bio/bio_omics/src/pipelines

# put IMGT data to database
python ${PIPELINE_DIR}/ab_imgt_genelist.py &
P1=$!
python ${PIPELINE_DIR}/ab_imgt_gene.py &
P2=$!
wait $P1 $P2

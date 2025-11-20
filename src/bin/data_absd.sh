#! /usr/bin/bash

PIPELINE_DIR=/home/yuan/bio/bio_omics/src/pipelines

echo "Try to put data of database ABSD into database"
parallel -j4 python ${PIPELINE_DIR}/absd_source.py ::: {0..4} &
P1=$!
parallel -j4 python ${PIPELINE_DIR}/absd_seq.py ::: {0..4} &
P2=$!
wait $P1 $P2

parallel -j8 python ${PIPELINE_DIR}/absd_faa.py ::: {0..360}

# igblastp
cd /home/yuan/data/ncbi-igblast-1.22.0 ; \
parallel -j8 python ${PIPELINE_DIR}/absd_igblastp.py ::: {0..360}

# then run the remaining pipelines
parallel -j8 python ${PIPELINE_DIR}/absd_align_imgt.py ::: {0..360}
parallel -j8 python ${PIPELINE_DIR}/absd_align_vregion.py ::: {0..360}


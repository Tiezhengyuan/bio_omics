#! /usr/bin/bash

PIPELINE_DIR=/home/yuan/bio/bio_omics/src/pipelines

echo "Try to download data from database AACDB"
python ${PIPELINE_DIR}/ab_aacdb.py

echo "Try to download data from database SAbDab"
python ${PIPELINE_DIR}/ab_sabdab.py
python ${PIPELINE_DIR}/ab_sabdab_insert.py

echo "Try to download data from database abYbank"
python ${PIPELINE_DIR}/ab_abybank.py

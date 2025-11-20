#! /usr/bin/bash

# NOTE: IMGT data should be available
#  bash src/bin/data_imgt.sh

PIPELINE_DIR=/home/yuan/bio/bio_omics/src/pipelines

# put IMGT data to database
# file genelist
python ${PIPELINE_DIR}/imgt_genelist.py &
P1=$!
# fasta files
python ${PIPELINE_DIR}/imgt_gene.py &
P2=$!
python ${PIPELINE_DIR}/imgt_geneseq.py &
P3=$!
# INN data
python ${PIPELINE_DIR}/imgt_inn.py &
P4=$!
wait $P1 $P2 $P3 $P4

# put chain firstly
python ${PIPELINE_DIR}/imgt_inn_chain.py
# Then add details by chain/region/domain
python ${PIPELINE_DIR}/imgt_inn_region.py &
P7=$!
python ${PIPELINE_DIR}/imgt_inn_cdomain.py &
P8=$!
python ${PIPELINE_DIR}/imgt_inn_vdomain.py &
P9=$!
python ${PIPELINE_DIR}/imgt_inn_cdr.py &
P10=$!
wait $P7 $P8 $P9 $P10

# identify INN sequences
python ${PIPELINE_DIR}/ab_blastp_db.py
python ${PIPELINE_DIR}/imgt_inn_blastp.py

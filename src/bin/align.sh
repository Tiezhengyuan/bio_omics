#! /usr/bin/bash

PIPELINE_DIR=/home/yuan/bio/bio_omics/src/pipelines

# align blastp
parallel -j4 python ${PIPELINE_DIR}/blastp_uniprot.py ::: {0..117}
parallel -j4 python ${PIPELINE_DIR}/align_uniprot.py ::: {0..117}

# align igblastp
cd /home/yuan/output/ncbi-igblast-1.22.0 ; \
parallel -j4 python ${PIPELINE_DIR}/igblastp_vregion.py ::: {0..117}
# then run the remaining pipelines
parallel -j4 python ${PIPELINE_DIR}/align_vregion.py ::: {0..117}
parallel -j4 python ${PIPELINE_DIR}/align_vfrag.py ::: {0..117}
# multiple alignment with IMGT germline
# parallel -j4 python ${PIPELINE_DIR}/imgt_msa.py ::: {0..117}


# identify INN sequences
# python ${PIPELINE_DIR}/ab_blastp_db.py
# python ${PIPELINE_DIR}/imgt_inn_blastp.py

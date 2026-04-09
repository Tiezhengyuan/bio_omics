#! /usr/bin/bash
APP=/home/yuan/bio/bio_omics/app.py

echo 'Data from IMGT/INN...'
# python ${APP} protein_inn
# identify INN sequences
# python ${APP} ab_blastp_db.py
# python ${APP} imgt_inn_blastp.py




# multiple alignment with IMGT germline
# parallel -j4 python ${APP} imgt_msa ::: {0..237}


echo "Try to put data of database ABSD into database"
parallel -j4 python ${PIPELINE_DIR}/absd_source.py ::: {0..4} &
P1=$!
parallel -j4 python ${PIPELINE_DIR}/absd_seq.py ::: {0..4} &
P2=$!
wait $P1 $P2

parallel -j8 python ${PIPELINE_DIR}/absd_faa.py ::: {0..360}


# igblastp
cd /home/yuan/data/ncbi-igblast-1.22.0 ; \
parallel -j8 python app.py absd_igblastp ::: {0..360}

# then run the remaining pipelines
parallel -j8 python app.py absd_align_imgt ::: {0..360}
parallel -j8 python app.py absd_align_vregion ::: {0..360}




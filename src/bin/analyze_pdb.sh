#! /usr/bin/bash

PIPELINE_DIR=/home/yuan/bio/bio_omics/src/pipelines

# align blastp
parallel -j4 python ${PIPELINE_DIR}/pdb_blastp_uniprot.py ::: {0..46}
parallel -j4 python ${PIPELINE_DIR}/put_align_uniprot.py ::: {0..46}

# align igblastp
cd /home/yuan/output/ncbi-igblast-1.22.0 ; \
parallel -j4 python ${PIPELINE_DIR}/pdb_igblastp.py ::: {0..46}
# then run the remaining pipelines
parallel -j4 python ${PIPELINE_DIR}/put_align_imgt.py ::: {0..46}
parallel -j4 python ${PIPELINE_DIR}/put_align_vregion.py ::: {0..46}
# multiple alignment with IMGT germline
parallel -j4 python ${PIPELINE_DIR}/imgt_msa.py ::: {0..46}


# atom distance
parallel -j4 python ${PIPELINE_DIR}/pdb_distance.py ::: {0..46}
# secondary structure
parallel -j4 python ${PIPELINE_DIR}/pdb_dssp.py ::: {0..46}
# require meta_combo_pdb
parallel -j4 python ${PIPELINE_DIR}/pdb_sasa.py ::: {0..46}

# identify INN sequences
python ${PIPELINE_DIR}/ab_blastp_db.py
python ${PIPELINE_DIR}/imgt_inn_blastp.py


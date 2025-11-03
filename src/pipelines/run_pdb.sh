#! /usr/bin/bash

PIPELINE_DIR=/home/yuan/bio/bio_omics/src/pipelines
LIST_DIR=/home/yuan/data/pdb/_list
n=46

# synchronize PDB
python ${PIPELINE_DIR}/pdb_rsync.py

# build table pdb
python ${PIPELINE_DIR}/put_pdb.py
python ${PIPELINE_DIR}/put_meta_pdb.py

# build pdb related tables
# pdb should be ready
parallel -j8 python ${PIPELINE_DIR}/put_chain.py ::: {0..${n}}

python ${PIPELINE_DIR}/put_compound.py &
P2=$!
python ${PIPELINE_DIR}/put_method.py &
P3=$!
python ${PIPELINE_DIR}/put_keyword.py &
P4=$!
wait $P2 $P3 $P4

# export chains in fasta or pdb
parallel -j4 python ${PIPELINE_DIR}/export_faa.py ::: {0..${n}}
parallel -j4 python ${PIPELINE_DIR}/export_chain_faa.py ::: {0..${n}}
parallel -j4 python ${PIPELINE_DIR}/export_chain_pdb.py ::: {0..${n}}

# calculate b-factor
parallel -j4 python ${PIPELINE_DIR}/pdb_bfactor.py ::: {0..${n}}
# chain pdb should be ready
parallel -j4 python ${PIPELINE_DIR}/pdb_dssp.py ::: {0..${n}}
# align blastp
parallel -j4 python ${PIPELINE_DIR}/pdb_blastp.py ::: {0..${n}}

# align igblastp
cd /home/yuan/data/ncbi-igblast-1.22.0
parallel -j4 python ${PIPELINE_DIR}/pdb_igblastp.py ::: {0..${n}}
parallel -j4 python ${PIPELINE_DIR}/put_align_imgt.py ::: {0..${n}}
parallel -j4 python ${PIPELINE_DIR}/put_align_vregion.py ::: {0..${n}}
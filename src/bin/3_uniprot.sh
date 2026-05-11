#! /usr/bin/bash
APP=/home/yuan/bio/bio_omics/app.py

# load uniprot_acc
python ${APP} uniprot

# load uniport annotations
python ${APP} uniprot_keyword &
P1=$!
python ${APP} uniprot_dbxrefs &
P2=$!
python ${APP} uniprot_feature &
P3=$!
wait $P1 $P2 $P3

echo 'load sequence to proseq_400, proseq_1000 or prseq_long'
python ${APP} uniprot_pro_seq


#alphafold prediction
python ${APP} uniprot_alphafold


# pLDDT of alphafold modelling
# python ${APP} uniprot_alphafold_plddt

# RMSD between chain pdb and uniprot prediction by alphafold
# parallel -j4 python ${APP} uniprot_alphafold_rmsd ::: {0..46}






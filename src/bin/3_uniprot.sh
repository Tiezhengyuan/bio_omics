#! /usr/bin/bash
$ bash src/bin/3_uniprot.sh
APP="uv --directory /home/yuan/bio/bio_omics run bomics"

# download data
${APP} download_uniprot

# load uniprot_acc
${APP} uniprot

# load uniport annotations
${APP} uniprot_keyword &
P1=$!
${APP} uniprot_dbxrefs &
P2=$!
${APP} uniprot_feature &
P3=$!
wait $P1 $P2 $P3

echo 'load sequence to proseq_400, proseq_1000 or prseq_long'
${APP} uniprot_pro_seq


#alphafold prediction
# ${APP} uniprot_af

# RMSD between chain pdb and uniprot prediction by alphafold
# parallel -j8 ${APP} uniprot_af_rmsd ::: {0..241}


# pLDDT of alphafold modelling
# python ${APP} uniprot_alphafold_plddt


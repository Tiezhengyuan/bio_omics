#! /usr/bin/bash
APP=/home/yuan/bio/bio_omics/app.py

echo "Align with reference of uniprot using blastp"
# pdb data, blastp, uniprot references are ready
parallel -j8 python ${APP} blastp_uniprot ::: {0..237}
parallel -j8 python ${APP} align_uniprot ::: {0..237}

echo "IgBlastp alignment on V-region and D/J region"
# PDB chain seq ready
# iglbastp+V-region refereces
# blastp + D/J region references
parallel -j8 python ${APP} meta_vdj ::: {0..237}

echo "Retrieve igblastp alignment results"
parallel -j8 python ${APP} align_vfrag ::: {0..237}
parallel -j8 python ${APP} align_vregion ::: {0..237}
parallel -j8 python ${APP} align_dregion ::: {0..237}
parallel -j8 python ${APP} align_jregion ::: {0..237}

parallel -j8 python ${APP} label_antibody ::: {0..237}


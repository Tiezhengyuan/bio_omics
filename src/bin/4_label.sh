#! /usr/bin/bash
APP=/home/yuan/bio/bio_omics/app.py

#### blastp

echo "Align against reference of uniprot using blastp"
# pdb data, blastp, uniprot references are ready
parallel -j8 python ${APP} meta_blastp5 ::: {0..240}
parallel -j8 python ${APP} align_uniprot ::: {0..240}

#### label antibody

echo "IgBlastp alignment against V-region and D/J region"
# PDB chain seq ready
# iglbastp+V-region refereces
# blastp + D/J region references
parallel -j8 python ${APP} meta_vdj ::: {0..240}
echo "Retrieve igblastp alignment results"
parallel -j8 python ${APP} align_vfrag ::: {0..240}
parallel -j8 python ${APP} align_vregion ::: {0..240}
parallel -j8 python ${APP} align_dregion ::: {0..240}
parallel -j8 python ${APP} align_jregion ::: {0..240}

# label chains
parallel -j8 python ${APP} label_antibody ::: {0..240}

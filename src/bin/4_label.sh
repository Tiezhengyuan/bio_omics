#! /usr/bin/bash
$ bash src/bin/4_label.sh
APP="uv --directory /home/yuan/bio/bio_omics run bomics"


##### step into ailgn 
#### blastp

echo "Align against reference of uniprot using blastp"
# pdb data, blastp, uniprot references are ready
parallel -j8 ${APP} meta_blastp5 ::: {0..241}
parallel -j8 ${APP} align_uniprot ::: {0..241}

#### label antibody

echo "IgBlastp alignment against V-region and D/J region"
# PDB chain seq ready
# iglbastp+V-region refereces
# blastp + D/J region references
parallel -j8 ${APP} meta_vdj ::: {0..241}
echo "Retrieve igblastp alignment results"
parallel -j8 ${APP} align_vfrag ::: {0..241}
parallel -j8 ${APP} align_vregion ::: {0..241}
parallel -j8 ${APP} align_dregion ::: {0..241}
parallel -j8 ${APP} align_jregion ::: {0..241}

# label chains
parallel -j8 ${APP} label_antibody ::: {0..241}



#! /usr/bin/bash
# bash src/bin/4_align.sh
APP="uv --directory /home/yuan/bio/bio_omics run bomics"


##### step into ailgn 
#### blastp

echo "Align against reference of uniprot using blastp"
# pdb data, blastp, uniprot references are ready
parallel -j8 ${APP} meta_blastp5 ::: {0..241}
parallel -j8 ${APP} align_uniprot ::: {0..241}

#### antibody region alignment: moved to antibody_omics (bin/6_align.sh)
# align_vdj, align_vfrag/vregion/dregion/jregion, label_antibody



#! /usr/bin/bash
APP="/home/yuan/bio/bio_omics/app.py"

echo "Push PDB to table pdb"
python ${APP} pdb
# python ${PIPELINE_DIR}/put_pdb2

echo "Prepare pdb_files into table meta_pdb"
python ${APP} meta_pdb
# python ${APP} put_meta_mmcif

# require meta_pdb
parallel -j8 python ${APP} pdb_compound ::: {0..237}
parallel -j8 python ${APP} pdb_method ::: {0..237}
parallel -j8 python ${APP} pdb_keyword ::: {0..237}

echo "Retriev B-factor..."
parallel -j8 python ${APP} meta_bfactor ::: {0..237}
parallel -j8 python ${APP} pdb_bfactor ::: {0..237}

echo "Analyze secondary structures"
parallel -j8 python ${APP} pdb_dssp ::: {0..237}

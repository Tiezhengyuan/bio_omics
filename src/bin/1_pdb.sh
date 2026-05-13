#! /usr/bin/bash
APP="/home/yuan/bio/bio_omics/app.py"

echo "Download PDB data to the local..."
python ${APP} download_pdb

echo "Push PDB to table pdb"
python ${APP} pdb
# python ${PIPELINE_DIR}/put_pdb2

echo "Detect pdb/cif file path into table meta_pdb"
python ${APP} meta_pdb

# require meta_pdb
parallel -j8 python ${APP} pdb_compound ::: {0..240}
parallel -j8 python ${APP} pdb_method ::: {0..240}
parallel -j8 python ${APP} pdb_keyword ::: {0..240}

#download original proseq in fasta
parallel -j8 python ${APP} pdb_protein ::: {0..240}

echo "Retriev B-factor..."
parallel -j8 python ${APP} meta_bfactor ::: {0..240}
parallel -j8 python ${APP} pdb_bfactor ::: {0..240}

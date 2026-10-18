#! /usr/bin/bash
set -e
APP="uv --directory /home/yuan/bio/bio_omics run bomics"

echo "Download PDB data to the local..."
${APP} download_pdb
${APP} download_cif

echo "Push PDB to table pdb"
${APP} pdb

echo "Detect pdb/cif file path into table meta_pdb"
${APP} meta_pdb

# require meta_pdb
parallel -j8 ${APP} pdb_compound ::: {0..241}
parallel -j8 ${APP} pdb_method ::: {0..241}
parallel -j8 ${APP} pdb_keyword ::: {0..241}

#download original proseq in fasta
parallel -j8 ${APP} pdb_protein ::: {0..241}

echo "Retriev B-factor..."
parallel -j8 ${APP} meta_bfactor ::: {0..241}
parallel -j8 ${APP} pdb_bfactor ::: {0..241}

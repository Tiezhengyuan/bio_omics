#! /usr/bin/bash
set -e
APP="uv --directory /home/yuan/bio/bio_omics run bomics"

echo "Download PDB data to the local..."
${APP} download_pdb
# ${APP} download_cif

echo "Push PDB to table pdb"
${APP} pdb

echo "Detect pdb/cif file path into table meta_pdb"
${APP} meta_pdb

# require meta_pdb
echo "Retrieve annotations of compounds, method, and keywords..."
parallel -j8 ${APP} pdb_compound ::: {0..258}
parallel -j8 ${APP} pdb_method ::: {0..258}
parallel -j8 ${APP} pdb_keyword ::: {0..258}

echo "Download original proseq in fasta from PDB..."
parallel -j8 ${APP} pdb_protein ::: {0..258}

echo "Retriev B-factor..."
parallel -j8 ${APP} meta_bfactor ::: {0..258}
parallel -j8 ${APP} pdb_bfactor ::: {0..258}

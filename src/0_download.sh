#! /usr/bin/bash
APP="/home/yuan/bio/bio_omics/app.py"

echo "Synchronize raw data from PDB"
python ${APP} download_pdb
# python ${APP} download_mmcif
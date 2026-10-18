#!/usr/bin/bash

output_pdb_dir="/home/yuan/omics/pdb"

# Check for any files (including hidden ones)

find /home/yuan/omics/output/pdb/combo5 -type d -empty -print -exec rm -d {} +

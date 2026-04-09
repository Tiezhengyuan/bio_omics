#! /usr/bin/bash
APP="/home/yuan/bio/bio_omics/app.py"

# identify chains
# parallel -j8 python ${APP} put_chain ::: {0..237}
parallel -j8 python ${APP} pdb_chainid ::: {0..237}

# retrieve sequences by chain
# require pdb_chain
parallel -j8 python ${APP} chain_seq ::: {0..237}
# update pdb_chain_seq
python ${APP} unique_chain_seq

#export sequences in faa
parallel -j8 python ${APP} pdb_faa ::: {0..237}
parallel -j8 python ${APP} meta_chain_faa ::: {0..237}

#identify unique sequences across all DBs
# insert data into protein_source
python ${APP} protein_pdb

# identify chain type
parallel -j8 python ${APP} pdb_chain_type ::: {0..237}

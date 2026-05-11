#! /usr/bin/bash
APP="/home/yuan/bio/bio_omics/app.py"

# identify chains
parallel -j8 python ${APP} pdb_chainid ::: {0..240}

# retrieve sequences by chain
# require pdb_chain
parallel -j8 python ${APP} chain_seq ::: {0..240}
# update pdb_chain_seq
python ${APP} unique_chain_seq

#export sequences in faa
parallel -j8 python ${APP} pdb_faa ::: {0..240}
parallel -j8 python ${APP} chain_faa ::: {0..240}

#parse pdb_chain_seq.chain_seq with pro_seq.seq_id
python ${APP} pdb_pro_seq

# identify chain type
parallel -j8 python ${APP} pdb_chain_type ::: {0..240}

# need pdb_chain_seq and pdb_chainid
parallel -j8 python ${APP} chain_pdb ::: {0..240}
# need chain_pdb
parallel -j8 python ${APP} chain_dssp ::: {0..240}

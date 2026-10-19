#! /usr/bin/bash
set -e
# bash bin/2_pdb_chain.sh
APP="uv --directory /home/yuan/bio/bio_omics run bomics"

# identify chains
parallel -j8 ${APP} pdb_chainid ::: {0..241}

# retrieve sequences by chain
parallel -j8 ${APP} pdb_chain_seq ::: {0..241}
${APP} unique_chain_seq

#export sequences in faa
parallel -j8 ${APP} pdb_faa ::: {0..241}

#parse pdb_chain_seq.chain_seq with pro_seq.seq_id
${APP} pdb_proseq

# identify chain type
parallel -j8 ${APP} pdb_chain_type ::: {0..241}

# need pdb_chain_seq and pdb_chainid
parallel -j8 ${APP} chain_pdb ::: {0..241}

#! /usr/bin/bash

# synchronize raw data from PDB
python app.py download_pdb
# python app.py download_mmcif

# PDB to pdb
python app.py pdb
# python ${PIPELINE_DIR}/put_pdb2

# pdb_files of PDB
python app.py meta_pdb
# python app.py put_meta_mmcif

# require meta_pdb
parallel -j8 python app.py pdb_compound ::: {0..237}
parallel -j8 python app.py pdb_method ::: {0..237}
parallel -j8 python app.py pdb_keyword ::: {0..237}
# parallel -j8 python app.py put_chain ::: {0..237}
parallel -j8 python app.py pdb_chainid ::: {0..237}

# require pdb_chain
parallel -j8 python app.py pdb_chain_seq ::: {0..237}
parallel -j8 python app.py unique_chain_seq
parallel -j8 python app.py pdb_faa ::: {0..237}
parallel -j8 python app.py meta_chain_faa ::: {0..237}
parallel -j8 python app.py meta_chain_pdb ::: {0..237}


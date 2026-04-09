

##############create and update
# 0. synchronize PDB
bash src/bin/0_download.sh

# 1. retrieve data
bash src/bin/1_pdb.sh

# 2. chain sequences
bash src/bin/2_chain.sh

# 3. annotate chain: alignment
bash src/bin/3_align.sh

# 4. split chains and relax
bash src/bin/4_simulate.sh

# 5. predict structures
bash src/bin/5_predict.sh


###############create only
#IMGT
bash src/bin/IMGT.sh
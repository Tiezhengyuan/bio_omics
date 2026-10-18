


#connect
alias docker-db="cd /home/yuan/bio/bio_omics; docker compose up -d"
#alias db="docker exec -it bomics-db mysql -u admin -p"
alias_db="title mysql;mysql -u admin -padmin COMPLEX2"

# synchronize coding
alias sync-bomics="rsync -av --exclude='__pycache__/' --exclude='.git/' --exclude='.venv/' yuan@am5:~/bio/bio_omics ~/bio/"
alias sync-abomics="rsync -av --exclude='__pycache__/' --exclude='.git/' --exclude='.venv/' yuan@am5:~/bio/antibody_omics ~/bio/"
alias sync-tools="rsync -av --exclude='__pycache__/' --exclude='.git/' yuan@am5:~/bio/tools ~/bio/"

# virtual env
alias env-bomics="cd ~/bio_omics; source .venv/bin/activate; title bomics"
alias env-abomics="cd ~/antibody_omics; source .venv/bin/activate; title abomics"
#alias mabs="conda activate simulate2;cd /home/yuan/bio/design_mabs"
#alias bind="conda activate bindcraft; cd /home/yuan/bio/binder/binder-design"

# run pipeline
# alias npr1="conda activate simulate; for i in {1..10}; do python ~/bio/design_mabs/NPR/design/npr1_predict.py; done; conda deactivate"
# alias npr2="conda activate simulate; for i in {1..10}; do python ~/bio/design_mabs/NPR/design/npr2_predict.py; done; conda deactivate"
# alias pred1="conda activate simulate2; cd ~/bio/bio_omics; python app.py chain_af"
# alias pred2="conda activate simulate2; cd ~/bio/bio_omics; python app.py combo2_af"
# alias pred2="conda activate simulate2; cd ~/bio/bio_omics; python app.py combo3_af"

#backup
alias bk-db="mysqldump COMPLEX2 > /media/yuan/data/bk_mysql/complex2_$(date +%Y%m%d).sql"
alias recover-db="docker exec -i bomics-db mysql -u admin -padmin COMPLEX2 < /media/yuan/data/bk_mysql/complex2_20260625.sql"
alias bk-bio="rsync -avhP -delete /home/yuan/bio/* /media/yuan/data/bk_bio/"
alias bk-predict="rsync -avhP -delete /home/yuan/omics/output/predict /media/yuan/data/bk_output/"

# upload results
alias up-predict="rsync -avzh -P --remove-source-files /home/yuan/output/predict/* yuan@i9:/home/yuan/output/predict/"
alias up-binder="rsync -avczh -P /home/yuan/output/predict/bindcarft yuan@i9:/home/yuan/output/predict/"

#clean
alias clean-gpu="killall -9 parallel || true && killall -9 python || true && killall -9 colabfold_batch || true;nvidia-smi"

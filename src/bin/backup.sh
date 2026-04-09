#!/usr/bin/bash

# example: 
# bash src/bin/backup.sh output-predict

name=$1
BK_DIR="/home/yuan/output/backup_mysql"
RSYNC_DIR="/media/yuan/complex_data"
OUTPUT_DIR="/home/yuan/output"
RAWDATA_DIR="/home/yuan/rawdata"

# backup mysql.complex2
if [ "$name" = "mysql" ]; then
    echo "Try to dump databased COMPLEX2 to ${BK_DIR}..."
    mysqldump -u admin COMPLEX2 > ${BK_DIR}/complex2_$(date +%Y%m%d).sql
    rsync -avzh --progress ${BK_DIR}/* ${RSYNC_DIR}/output/mysql/

# synchronize mysql backup to portable hard drive
elif [ "$name" = "output-pdb" ]; then
    echo "Try to ${OUTPUT_DIR}/pdb to ${RSYNC_DIR}..."
    rsync -avzh --progress --delete ${OUTPUT_DIR}/pdb ${RSYNC_DIR}/output/
elif [ "$name" = "output-predict" ]; then
    echo "Try to backup ${OUTPUT_DIR}/predict to ${RSYNC_DIR}..."
    rsync -avzh --progress --delete ${OUTPUT_DIR}/predict ${RSYNC_DIR}/output/
else
    echo "arguments should be mysql|rsync-mysql|pdb-data."
fi

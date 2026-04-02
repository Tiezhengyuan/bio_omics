#!/usr/bin/bash

# example: bash src/bin/backup.sh <data>

name=$1
BK_DIR="/home/yuan/output/backup_mysql"
RSYNC_DIR="/media/yuan/complex"
DATA_DIR="/home/yuan/output/pdb"

# backup mysql.complex2
if [ "$name" = "mysql" ]; then
    echo "Try to dump databased COMPLEX2 to ${BK_DIR}..."
    mysqldump -u admin COMPLEX2 > ${BK_DIR}/complex2_$(date +%Y%m%d).sql
# synchronize mysql backup to portable hard drive
elif [ "$name" = "rsync-mysql" ]; then
    echo "Try to dump databased COMPLEX2 to ${RSYNC_DIR}..."
    rsync -avz --progress ${BK_DIR} ${RSYNC_DIR}/
elif [ "$name" = "pdb-data" ]; then
    echo "Try to backup output/pdb to ${RSYNC_DIR} in portable hard drive..."
    rsync -avz --progress ${DATA_DIR} ${RSYNC_DIR}/pdb_data
else
    echo "arguments should be mysql|rsync-mysql|pdb-data."
fi

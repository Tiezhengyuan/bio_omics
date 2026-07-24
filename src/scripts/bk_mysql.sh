#!/usr/bin/bash

BK_DIR="/home/yuan/output/backup_mysql"
RSYNC_DIR="/media/yuan/complex_data"

echo "Try to dump databased COMPLEX2 to ${BK_DIR}..."
mysqldump -u admin COMPLEX2 > ${BK_DIR}/complex2_$(date +%Y%m%d).sql
rsync -avzh --progress ${BK_DIR}/* ${RSYNC_DIR}/output/mysql/

#!/usr/bin/bash

# example: bash src/bin/backup_mysql.sh

BK_DIR="/home/yuan/output/backup_mysql"

mysqldump -u admin -p COMPLEX2 > ${BK_DIR}/complex2_$(date +%Y%m%d).sql
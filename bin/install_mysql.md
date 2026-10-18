
# native mysql
```
sudo apt --update
sudo apt  --upgrade
```

```
sudo apt install mysql-server
sudo systemctl start mysql
sudo systemctl enable msyql
```

secure installation

```
sudo mysql_secure_installation
```

login
```
sudo mysql
```

```
SET GLOBAL validate_password.policy = LOW;
SET GLOBAL validate_password.length = 4;

-- 1. Create the local access user
CREATE USER 'admin'@'localhost' IDENTIFIED BY 'admin';

-- 2. Create the remote access user (% allows any IP, or replace % with a specific IP)
CREATE USER 'admin'@'%' IDENTIFIED BY 'admin';

-- 3. Grant privileges to the local user (replace 'mydatabase' with your actual DB name)
GRANT ALL PRIVILEGES ON *.* TO 'admin'@'localhost';

-- 4. Grant privileges to the remote user
GRANT ALL PRIVILEGES ON *.* TO 'admin'@'%';

-- 5. Reload the grant tables to apply changes
FLUSH PRIVILEGES;
```

```
mysql -u admin -padmin COMPLEX2
```

check port
``
sudo ss -tuln | grep 3306
``

socket, ~/.my.cnf

[client]
socket=/var/run/mysqld/mysqld.sock



/etc/mysql/mysql.conf.d/mysqld.cnf

bind-address = 0.0.0.0

# deploy msyql data by docker

```
cd ~/bio/bio_omics
docker compose up -d
```

login database
```
docker exec -it bomics-db mysql -u admin -padmin COMPLEX2
```

stop container
```
docker compose down -v
```
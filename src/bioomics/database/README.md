


## install on Ubuntu Linux

```
sudo apt-get update
sudo apt-get install mysql-server
systemctl status mysql
```

configuration
```
sudo mysql_secure_installation
```

```
sudo mysql
# sudo mysql -uroot

```


```
select user, host, plugin from mysql.user;
```

change 
```
ALTER USER 'root'@'localhost' IDENTIFIED WITH mysql_native_password BY 'root';
FLUSH PRIVILEGES;
```

log in database
```
mysql -u root -p
```

```
mysql --version
```

```
mysql -u root -p protein < iedb_public.sql
```



uninstall mysql
```
sudo apt-get remove mysql-client mysql-server -y
sudo apt-get autoremove -y
sudo apt-get autoclean -y

#remove residuals
#sudo mv /var/lib/mysql /var/lib/mysql_directory_backup
sudo rm -fr /var/lib/mysql

# remove users
sudo deluser --remove-home mysql
sudo delgroup mysql
```
import mysql.connector

class BuildDb:

    def __init__(self, data:list=None, verbose:bool=False, chunk_size:int=None):
        self.data = data
        self.verbose = verbose
        self.chunk_size = int(chunk_size) if chunk_size else 100
        self.config = {}

    def execute_sql(self, query):
        try:
            connection = mysql.connector.connect(
                host=self.config['host'],
                user=self.config['user'],
                password=self.config['password'],
            )
            cursor = connection.cursor(dictionary=True)
            cursor.execute(query)
            connection.commit()
            return cursor
        except Exception as e:
            if self.verbose:
                print(f"Error: {e}")
        finally:
            if 'connection' in locals() and connection.is_connected():
                cursor.close()
                connection.close()
        return None

    def drop_database(self, db_name=None):
        db_name = db_name if db_name else self.config['database']
        query = f"DROP DATABASE {db_name};"
        result = self.execute_sql(query)
        if self.verbose:
            if result:
                print(f"Database {db_name} is dropped.")

    def create_database(self, db_name=None):
        db_name = db_name if db_name else self.config['database']
        query = f"CREATE DATABASE {db_name};"
        result = self.execute_sql(query)
        if self.verbose:
            if result:
                print(f"Database {db_name} is created.")

    def execute_db_sql(self, query):
        try:
            connection = mysql.connector.connect(
                host=self.config['host'],
                user=self.config['user'],
                password=self.config['password'],
                database=self.config['database'],
            )
            cursor = connection.cursor(dictionary=True)
            cursor.execute(query)
            connection.commit()
            return cursor
        except Exception as e:
            if self.verbose:
                print(f"Error: {e}")
        finally:
            if 'connection' in locals() and connection.is_connected():
                cursor.close()
                connection.close()
        return None

    def drop_table(self, table_name):
        query = f"DROP TABLE {table_name};"
        return self.execute_db_sql(query)

    def create_table(self, query):
        try:
            connection = mysql.connector.connect(**self.config)
            cursor = connection.cursor(dictionary=True)
            cursor.execute(query)
            connection.commit()
            return cursor
        except Exception as e:
            if self.verbose:
                print(f"Error: {e}")
        finally:
            if 'connection' in locals() and connection.is_connected():
                cursor.close()
                connection.close()
        return None

    def empty_table(self, table_name):
        if self.verbose:
            print(f"Try to delete all records from table {table_name}...")
        query = f"DELETE FROM {table_name};"
        return self.execute_db_sql(query)

    def insert_records(self, insert_query, records):
        try:
            connection = mysql.connector.connect(**self.config)
            cursor = connection.cursor(dictionary=True)
            cursor.executemany(insert_query, records)
            connection.commit()
            return cursor.rowcount
        except Exception as e:
            if self.verbose:
                print(f"Error: {e}. data={records}")
            pass
        finally:
            if 'connection' in locals() and connection.is_connected():
                cursor.close()
                connection.close()
        return None

    def insert_data(self, insert_query, data:list):
        n, m = 0, 0
        for rec in data:
            records = [tuple(rec),]
            result = self.insert_records(insert_query, records)
            if result:
                n+= 1
            else:
                m += 1
        print(f"table complex: total={len(data)}, succeed={n}, failed={m}")

    def insert(self, table_name, insert_query, func, pool_size:int=None, **kwargs):
        '''
        insert data into table
        '''
        pool_size = self.chunk_size if pool_size is None else int(pool_size)
        print(f"\nBegin to insert records into {table_name}. Scan Data...", end='')
        # get records: self.data is pdb paths
        n, succeed = 0, 0
        pool, recycle_pool = [], []
        for arg in self.data:
            n += 1
            if n % 1_000 == 0:
                print(n, end='..')
            # collect record: arg is pdb_file or dat_record
            records = func(arg, **kwargs)
            pool.extend(records)
            # insert records into table
            if len(pool) >= pool_size:
                m = self.insert_records(insert_query, pool)
                if m:
                    succeed += len(pool)
                else:
                    recycle_pool += pool
                pool, recycle_pool = [], []
        if pool:
            m = self.insert_records(insert_query, pool)
            if m:
                succeed += len(pool)
            else:
                recycle_pool += pool
        print(f"\nTry again: insert one by one...")
        failed = []
        for rec in recycle_pool:
            m = self.insert_records(insert_query, [rec,])
            if m:
                succeed += 1
            else:
                failed.append(rec)
        print(f"failed:{len(failed)}, succeed={succeed}.\n")
        return failed


'''
configurations for mysql.connector.connection
'''

import os
import mysql.connector

class Db:

    def __init__(self, verbose:bool):
        self.verbose = verbose
        self.config = {
            'host': os.getenv('MYSQL_HOST'),
            'user': os.getenv('MYSQL_USER'),
            'password': os.getenv('MYSQL_PW'),
            'use_pure': True,
            'buffered': True,
            'autocommit': True,    # avoid transaction issues
            'connection_timeout': 60,
            'pool_size': 1,        # Single connection for stability
        }

    def execute(self, query) -> bool:
        try:
            connection = mysql.connector.connect(**self.config)
            cursor = connection.cursor(dictionary=True)
            cursor.execute(query)
            connection.commit()
            return True
        except mysql.connector.Error as e:
            # Roll back if there is an error
            connection.rollback()
            if self.verbose:
                print(f"Error: {e}")
        except Exception as e:
            if self.verbose:
                print(f"Error: {e}")
        finally:
            if 'connection' in locals() and connection.is_connected():
                cursor.close()
                connection.close()
        return False

    def execute_fetch(self, query):
        try:
            connection = mysql.connector.connect(**self.config)
            cursor = connection.cursor(dictionary=True)
            cursor.execute(query)
            res = cursor.fetchall()
            connection.commit()
            return res
        except Exception as e:
            if self.verbose:
                print(f"Error: {e}")
        finally:
            if 'connection' in locals() and connection.is_connected():
                cursor.close()
                connection.close()
        return -1

    def execute_many(self, query, records) -> int:
        try:
            connection = mysql.connector.connect(**self.config)
            cursor = connection.cursor(dictionary=True)
            cursor.executemany(query, records)
            connection.commit()
            return cursor.rowcount
        except Exception as e:
            if self.verbose:
                print(f"Error: {e}")
        finally:
            if 'connection' in locals() and connection.is_connected():
                cursor.close()
                connection.close()
        return 0
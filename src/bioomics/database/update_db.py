from copy import deepcopy
import mysql.connector

from .db import Db

class UpdateDb(Db):

    def __init__(self, verbose:bool=True):
        super().__init__(verbose)

    def execute_update(self, query, record):
        '''
        execute one query at a time
        '''
        record = deepcopy(record)
        try:
            connection = mysql.connector.connect(**self.config)
            cursor = connection.cursor(dictionary=True)
            cursor.execute(query, record)
            connection.commit()
            return cursor.rowcount
        except Exception as e:
            if self.verbose:
                print(f"Error: {e}")
        finally:
            if 'connection' in locals() and connection.is_connected():
                cursor.close()
                connection.close()
        return None


import mysql.connector
from copy import deepcopy

class UpdateDb:

    def __init__(self, verbose:bool=False):
        self.verbose = verbose
        self.config = {}

    def execute_update(self, query, record):
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

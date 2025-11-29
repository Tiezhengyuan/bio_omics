from copy import deepcopy
import mysql.connector

from .db import Db

class UpdateDb(Db):

    def __init__(self, verbose:bool=False):
        super().__init__(verbose)

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

    def delete_rows(self, table_name, col_name, value):
        query = f"""
            DELETE FROM {table_name}
            WHERE {col_name} = %s
        """
        adr = (value,)
        return self.execute_update(query, adr)
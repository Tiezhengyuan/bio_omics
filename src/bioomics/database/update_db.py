from copy import deepcopy
import mysql.connector

from .db import Db

class UpdateDb(Db):

    def __init__(self, table_name, verbose:bool=False):
        self.table_name = table_name
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

    def delete_rows(self, col_name, values):
        placeholders = ', '.join(['%s',] * len(values))
        query = f"""
            DELETE FROM {self.table_name}
            WHERE {col_name} IN ({placeholders})
        """
        return self.execute_update(query, values)

    def delete_null_rows(self, col_name):
        query = f"""
            DELETE FROM {self.table_name}
            WHERE {col_name} IS NULL
        """
        return self.execute_update(query, None)
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

    # def empty_tables(self, table_name:str):
    #     query = f"DELETE FROM {table_name};"
    #     return self.execute_update(query, None)

    def delete_rows(self, table_name, col_name:str, values:list):
        if self.verbose:
            print(f"Try to delete rows ={values} from {table_name}-{col_name}")
        placeholders = ', '.join(['%s',] * len(values))
        query = f"""
            DELETE FROM {table_name}
            WHERE {col_name} IN ({placeholders})
        """
        return self.execute_update(query, values)

    def delete_null_rows(self, table_name, col_name):
        query = f"""
            DELETE FROM {table_name}
            WHERE {col_name} IS NULL
        """
        return self.execute_update(query, None)
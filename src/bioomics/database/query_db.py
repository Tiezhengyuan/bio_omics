import gc
import pandas as pd
import mysql.connector
import time

from .db import Db

class QueryDb(Db):

    def __init__(self, verbose:bool=False):
        super().__init__(verbose)

    def list_data(self, query:str, is_df:bool=False) -> list:
        for _ in range(5):
            try:
                with mysql.connector.connect(**self.config) as connection:
                    with connection.cursor(dictionary=True) as cursor:
                        cursor.execute(query)
                        rows = cursor.fetchall()
                        if is_df:
                            df = pd.DataFrame(rows)
                            return df
                        return rows
            except Exception as e:
                if self.verbose:
                    print(f"error: can't connect mysql. error={e}")
                gc.collect()
                time.sleep(5)
        return []

    # TODO: old version, deprecated in the future.
    def list_data0(self, query:str, is_df:bool=False):
        try:
            connection = mysql.connector.connect(**self.config)
            cursor = connection.cursor(dictionary=True)
            cursor.execute(query)
            rows = cursor.fetchall()
            if is_df:
                df = pd.DataFrame(rows)
                return df
            return rows
        except Exception as e:
            if self.verbose:
                print(f"Error: {e}")
        finally:
            if 'connection' in locals() and connection.is_connected():
                cursor.close()
                connection.close()
        return []

    def table_names(self):
        query = "SHOW TABLES;"
        res = self.list_data(query, False)
        res = pd.DataFrame(res).iloc[:,0]
        return tuple(res)

    def unique_values(self, table_name, col_name:str) -> list:
        query = f"""
            SELECT DISTINCT {col_name}
            FROM {table_name}
            WHERE {col_name} IS NOT NULL;
        """
        rows = self.list_data(query, True)
        if self.verbose:
            print(f"Detect {len(rows)} records in {table_name}.{col_name}.")
        if len(rows) > 0:
            return list(rows[col_name])
        return []

    def next_id(self, table_name, col_name):
        query = f"""
            SELECT MAX({col_name}) AS id FROM {table_name};
        """
        res = self.list_data(query, False)
        curr = res[0]['id']
        return int(curr) + 1 if curr else 1

    def max_value(self, table_name:str, col:str):
        query = f"""
            SELECT MAX({col}) AS val
            FROM {table_name}
        ;"""
        res = self.list_data(query, False)
        return res[0]['val']
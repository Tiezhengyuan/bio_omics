import pandas as pd
import mysql.connector

class QueryDb:

    def __init__(self, verbose:bool=False):
        self.verbose = verbose
        self.config = {}

    def list_data(self, query:str, is_df:bool=False):
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

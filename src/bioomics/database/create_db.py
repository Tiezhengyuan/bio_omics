'''
'''
from .db import Db

class CreateDb(Db):

    def __init__(self, verbose:bool=False):
        super().__init__(verbose)
    
    def create_table(self, table_name:str, cols_str:str):
        query = f"""
            CREATE TABLE {table_name} ({cols_str})
        """
        return self.execute(query)

    def create_database(self, db_name=None):
        db_name = db_name if db_name else self.config['database']
        query = f"CREATE DATABASE {db_name};"
        result = self.execute(query)
        if self.verbose:
            if result:
                print(f"Database {db_name} is created.")

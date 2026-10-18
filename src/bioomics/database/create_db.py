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
        db_name = db_name if db_name else self.config.get('database')
        if not db_name:
            return False
        query = f"CREATE DATABASE {db_name};"
        # the database doesn't exist yet: connect without selecting it
        database = self.config.pop('database', None)
        try:
            result = self.execute(query)
        finally:
            if database:
                self.config['database'] = database
        if self.verbose:
            if result:
                print(f"Database {db_name} is created.")
        return result

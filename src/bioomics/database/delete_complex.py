from collections import defaultdict
import os
import pandas as pd

from .delete_db import DeleteDb

class DeleteComplex(DeleteDb):

    def __init__(self, verbose:bool=None):
        super().__init__(verbose)
        self.config['database'] = os.getenv('MYSQL_DB', 'COMPLEX2')
    
    def remove_null_rows(self, table_name:str, col_name:str, del_pool:list):
        n = len(del_pool)
        num_del = 0
        if n > 0:
            print(f"Number of {n} rows are invalid.", end=' ')
            chunk_size = 100
            end = n - chunk_size if n > chunk_size else n
            for i in range(0, end, chunk_size):
                values = del_pool[i:i+chunk_size]
                res = self.delete_rows(table_name, col_name, values)
                if res:
                    num_del += len(values)
                    print(i, end=' ')
            print(f'{n}\n')
        else:
            print(f"PASS: Validate {table_name}.{col_name}.")
        return num_del
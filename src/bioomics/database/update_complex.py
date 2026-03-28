import os
from datetime import datetime

from .update_db import UpdateDb

class UpdateComplex(UpdateDb):

    def __init__(self, table_name:str, verbose:bool=False):
        super().__init__(table_name, verbose)
        self.config['database'] = os.getenv('MYSQL_DB', 'COMPLEX2')
    
    def pdb_path_null(self, data:dict, col_name=None):
        '''
        set relative_path NULL if a file doesn't exist
        '''
        if col_name is None:
            col_name = 'relative_path'
        update_query = f"""
            UPDATE {self.table_name}
            SET {col_name} = NULL
            WHERE pdb_id = %(pdb_id)s
        """
        if self.verbose:
            print(f"Try to update table: {data}...")
        result = self.execute_update(update_query, data)
        return result

    def chain_path_null(self, data:dict, col_name):
        '''
        set relative_path NULL if a file doesn't exist
        '''
        update_query = f"""
            UPDATE {self.table_name}
            SET {col_name} = NULL
            WHERE pdb_id = %(pdb_id)s AND chain_id = %(chain_id)s
        """
        if self.verbose:
            print(f"Try to update table: {data}...")
        result = self.execute_update(update_query, data)
        return result
    
    def combo_path_null(self, data:dict, col_name):
        '''
        set relative_path NULL if a file doesn't exist
        '''
        update_query = f"""
            UPDATE {self.table_name}
            SET {col_name} = NULL
            WHERE pdb_id = %(pdb_id)s AND combo_id = %(combo_id)s
        """
        if self.verbose:
            print(f"Try to update table: {data}...")
        result = self.execute_update(update_query, data)
        return result
    
    def update_compound(self, df, col_name:str):
        '''
        table is compound
        '''
        update_query = f"""
            UPDATE {self.table_name}
            SET {col_name} = %({col_name})s
            WHERE pdb_id = %(pdb_id)s AND compound_no = %(compound_no)s
        """
        n,m=0,0
        cols = [col_name, 'pdb_id', 'compound_no']
        data = df[cols].to_dict(orient='records')
        for rec in data:
            r = self.execute_update(update_query, rec)
            if r:
                n+= 1
            else:
                m += 1
        print(f"{col_name}: succeed={n}, unchanged={m}")
    
    def update_compound_table(self, col_name, df):
        '''
        The compound tables could be compountd, abag
        where primary key is )pdb_id, compound_no)
        '''
        update_query = f"""
            UPDATE {self.table_name}
            SET {col_name} = %({col_name})s
            WHERE pdb_id = %(pdb_id)s AND compound_no = %(compound_no)s
        """
        n,m=0,0
        cols = ['pdb_id', 'compound_no', col_name]
        data = df[cols].to_dict(orient='records')
        for rec in data:
            # records = [rec,]
            r = self.execute_update(update_query, rec)
            if r:
                n+= 1
            else:
                m += 1
        print(f"{self.table_name}.{col_name}: succeed={n}, unchanged={m}")

    def update_chain_seq(self, data:list):
        update_query = f"""
            UPDATE {self.table_name}
            SET first_chain_id = %s
            WHERE pdb_id = %s AND chain_id = %s
        """
        n,m=0,0
        for rec in data:
            r = self.execute_update(update_query, rec)
            if r:
                n+= 1
            else:
                m += 1
        print(f"{self.table_name}.first_chain_id: succeed={n}, unchanged={m}")
    
    def update_execute_time(self, args):
        t1 = datetime.now()
        query = f"""
            update {self.table_name}
            set end_at = '{t1}'
            where table_name = '{args['table_name']}'
                and pdb_group = {args['pdb_group']}
                and host = '{args['host']}'
        ;"""
        self.execute_update(query, None)
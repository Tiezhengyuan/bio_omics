from .update_db import UpdateDb
from .build_db import BuildDb
from .build_complex import BuildComplex
from .query_complex import QueryComplex

class UpdateComplex(BuildDb):

    def __init__(self, data, verbose:bool=False):
        super().__init__(data, verbose)
        self.config = {
            'host': 'localhost',
            'user': 'admin',
            'password': 'strong_password',
            'database': 'COMPLEX2'
        }

    def update_compound(self, df, col_name:str):
        '''
        table is compound
        '''
        update_query = f"""
            UPDATE compound
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
    
    def update_compound_table(self, table_name, col_name, df):
        '''
        The compound tables could be compountd, abag
        where primary key is )pdb_id, compound_no)
        '''
        update_query = f"""
            UPDATE {table_name}
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
        print(f"{table_name}.{col_name}: succeed={n}, unchanged={m}")


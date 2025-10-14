import os
import pandas as pd
from copy import deepcopy

from .query_db import QueryDb

class QueryComplex(QueryDb):

    def __init__(self, verbose:bool=None):
        super().__init__(verbose)
        self.config = {
            'host': 'localhost',
            'user': 'admin',
            'password': 'strong_password',
            'database': 'COMPLEX2'
        }

    def get_latest_date(self, table_name):
        query = f"""
            SELECT MAX(release_date) AS date FROM {table_name};
        """
        res = self.list_data(query, False)
        return res[0]['date']

    def get_unique_values(self, table_name, col_name:str) -> list:
        query = f"""
            SELECT DISTINCT {col_name}
            FROM {table_name}
            WHERE {col_name} is NOT NULL;
        """
        rows = self.list_data(query, True)
        if self.verbose:
            print(f"Detect {len(rows)} records in {table_name}.{col_name}.")
        if len(rows) > 0:
            return list(rows[col_name])
        return []

    def delete_pdb_ids(self, table_name:str):
        query = f"""
            DELETE FROM {table_name}
            WHERE chain_type = "x";
        """
        pdb_ids = self.list_data(query)
        return tuple(pdb_ids['pdb_id'])

    def count_complex(self):
        query = f"""
            SELECT complex_type, COUNT(DISTINCT pdb_id) AS pdb_ids
            FROM complex
            GROUP BY complex_type
            ORDER BY complex_type;
        """
        res = self.list_data(query)
        return pd.DataFrame(res)

    def scan_abag(self):
        query = f"""
            SELECT A.pdb_id, A.compound_no, A.chain_type, B.chain, C.local_path
            FROM abag A, compound B, pdb C
            WHERE A.pdb_id = B.pdb_id
                AND A.compound_no = B.compound_no
                AND A.pdb_id = C.pdb_id
            ORDER BY A.pdb_id, A.compound_no;
        """
        res = self.list_data(query)
        g = pd.DataFrame(res).groupby(['pdb_id', 'local_path'])
        for (pdb_id, local_path), sub in g:
            rec = {
                'pdb_id': pdb_id,
                'pdb_file': local_path,
                'antibody': [],
                'antigen': [],
            }
            for item in sub.to_dict(orient='records'):
                chain_type = item['chain_type']
                for chain_id in item['chain'].split(','):
                    _rec = (chain_type, chain_id.upper())
                    if chain_type == 'antigen':
                        rec['antigen'].append(_rec)
                    else:
                        rec['antibody'].append(_rec)
            yield rec

    def get_seq(self):
        query = """
            select A.pdb_id, A.compound_no, A.chain_id, A.chain_type, B.pro_seq
            from abag A, chain B
            where A.pdb_id = B.pdb_id
                AND A.compound_no = B.compound_no
                AND A.chain_id = B.chain_id
                AND B.pro_seq is NOT NULL
                AND A.chain_type IN ('light','heavy','single');
        """
        rows = self.list_data(query)
        print('sequencines of abag:', len(rows))
        res = {}
        for rec in rows:
            key = '|'.join([rec['pdb_id'], rec['chain_id'], rec['chain_type'],
                str(len(rec['pro_seq'])),])
            res[key] = str(rec['pro_seq'])
        return res

    
    def dup_uniprot_align(self):
        query = """
            WITH dup AS (
                SELECT pdb_id, model_id, chain_id, uniprot_acc
                FROM uniprot_align
                GROUP BY pdb_id, model_id, chain_id, uniprot_acc
                HAVING COUNT(*) >1)
            SELECT A.pdb_id, A.chain_id, A.uniprot_acc
            FROM uniprot_align A
            INNER JOIN dup B
            ON A.pdb_id=B.pdb_id AND A.model_id = B.model_id
             AND A.chain_id = B.chain_id
             AND A.uniprot_acc = B.uniprot_acc;
        """
        rows = self.list_data(query)
        print('duplicates of uniprot_align:', len(rows))
        return rows

    def antibody(self):
        query = """
            SELECT * FROM view_pdb_uniprot
            WHERE keyword = 'Immunoglobulin' OR keyword LIKE '%antibody%';
        """
        rows = self.list_data(query)
        return rows

    def antibody_complex(self):
        query = """
            SELECT * FROM view_pdb_uniprot
            WHERE pdb_id IN(
                SELECT DISTINCT pdb_id
                FROM view_pdb_uniprot
                WHERE keyword = 'Immunoglobulin' OR keyword LIKE '%antibody%'
            );
        """
        rows = self.list_data(query)
        return rows


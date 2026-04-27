import os
import pandas as pd

from .query_db import QueryDb

class QueryAntibody(QueryDb):
    def __init__(self, verbose:bool=False):
        super().__init__(verbose)
        self.config['database'] = os.getenv('MYSQL_DB', 'COMPLEX2')

    def map_chains(self) -> dict:
        query = """
            SELECT pdb_id, chain_no, chain_type
            FROM view_antibody
            WHERE model_no = 0
        ;"""
        res = {}
        for row in self.list_data(query):
            pdb_id = row['pdb_id']
            if pdb_id not in res:
                res[pdb_id] = {}
            chain_no, chain_type = row['chain_no'], row['chain_type']
            res[pdb_id][chain_no] = chain_type
        return res

    def antibody_chains(self, pdb_group:int, groupby:list=None):
        query = f"""
            SELECT pdb_id, model_no, chain_id, first_chain_id, chain_no,
                chain_type, relative_faa, chain_pdb
            FROM view_antibody
            WHERE pdb_group = {pdb_group}
                AND relative_faa IS NOT NULL
                AND chain_pdb IS NOT NULL
        ;"""
        df = self.list_data(query, True)
        if groupby:
            return df.groupby(groupby)
        return df

    def heavy_chains(self, pdb_group:int):
        df = self.antibody_chains(pdb_group)
        df = df[df['chain_type']=='H']
        return df
    
    def light_chains(self, pdb_group:int):
        df = self.antibody_chains(pdb_group)
        df = df[df['chain_type'].isin(['K', 'L'])]
        return df
    
    def target_chains(self, pdb_group:int):
        query = f"""
            SELECT pdb_id, model_no, chain_id, chain_no, relative_faa, chain_pdb
            FROM view_pdb_chain
            WHERE pdb_group = {pdb_group}
                AND relative_faa IS NOT NULL
                AND chain_pdb IS NOT NULL
                AND chain_id NOT IN (
                    SELECT chain_id FROM view_antibody
                )
        ;"""
        df = self.list_data(query, True)
        return df
    
    def group_to_dict(self, df):
        res = {}
        g = df.groupby(['pdb_id', 'model_no'])
        for key, sub in g:
            res[key] = sub
        return res
    
    def dimer_heavy_light(self, pdb_group:int) -> list:
        res = []
        g = self.antibody_chains(pdb_group, ['pdb_id', 'model_no'])
        for (pdb_id, model_no), sub in g:
            heavy = sub[sub['chain_type']=='H']
            light = sub[sub['chain_type'].isin(['K', 'L'])]
            if len(heavy) > 0 and len(light) > 0:
                for _, row1 in heavy.iterrows():
                    for _, row2 in light.iterrows():
                        chain_combo = [row1['chain_no'], row2['chain_no']]
                        chain_label = [row1['chain_type'], row2['chain_type']]
                        chain_ids = [row1['chain_id'], row2['chain_id']]
                        res.append({
                            'pdb_id': pdb_id,
                            'model_no': model_no,
                            'chain_combo': chain_combo,
                            'chain_label': chain_label,
                            'chain_ids': chain_ids,
                        })
        return res
    
    def dimer_antibody_target(self, pdb_group:int) -> list:
        res = []
        g1 = self.antibody_chains(pdb_group, ['pdb_id', 'model_no'])
        df2 = self.target_chains(pdb_group)
        dict2 = self.group_to_dict(df2)
        for key, sub1 in g1:
            if key in dict2:
                sub2 = dict2[key]
                for _, row1 in sub1.iterrows():
                    for _, row2 in sub2.iterrows():
                        chain_combo = [row1['chain_no'], row2['chain_no']]
                        chain_label = [row1['chain_type'], 'target']
                        chain_ids = [row1['chain_id'], row2['chain_id']]
                        res.append({
                            'pdb_id': key[0],
                            'model_no': key[1],
                            'chain_combo': chain_combo,
                            'chain_label': chain_label,
                            'chain_ids': chain_ids,
                        })
        return res

    def trimer_heavy_light_target(self, pdb_group:int) -> list:
        res = []
        df1 = self.heavy_chains(pdb_group)
        g1 = df1.groupby(['pdb_id', 'model_no'])
        df2 = self.light_chains(pdb_group)
        dict2 = self.group_to_dict(df2)
        df3 = self.target_chains(pdb_group)
        dict3 = self.group_to_dict(df3)
        for key, sub1 in g1:
            if key in dict2 and key in dict3:
                sub2, sub3 = dict2[key], dict3[key]
                for _, row1 in sub1.iterrows():
                    for _, row2 in sub2.iterrows():
                        for _, row3 in sub3.iterrows():
                            chain_combo = [row1['chain_no'], row2['chain_no'], row3['chain_no']]
                            chain_label = [row1['chain_type'], row2['chain_type'], 'target']
                            chain_ids = [row1['chain_id'], row2['chain_id'], row3['chain_id']]
                            res.append({
                                'pdb_id': key[0],
                                'model_no': key[1],
                                'chain_combo': chain_combo,
                                'chain_label': chain_label,
                                'chain_ids': chain_ids,
                            })
        return res

    def tetramer_antibody(self, pdb_group:int) -> list:
        res = []
        df1 = self.heavy_chains(pdb_group)
        g1 = df1.groupby(['pdb_id', 'model_no'])
        df2 = self.light_chains(pdb_group)
        dict2 = self.group_to_dict(df2)
        for key, sub1 in g1:
            if key in dict2:
                sub2 = dict2[key]
                if len(sub1) == 2 and len(sub2) == 2:
                    chain_combo = pd.concat([sub1['chain_no'], sub2['chain_no']])
                    chain_label = pd.concat([sub1['chain_type'], sub2['chain_type']])
                    chain_ids = pd.concat([sub1['chain_id'], sub2['chain_id']])
                    res.append({
                        'pdb_id': key[0],
                        'model_no': key[1],
                        'chain_combo': chain_combo,
                        'chain_label': chain_label,
                        'chain_ids': chain_ids,
                    })
        return res
    
    def pentamer_antibody_target(self, pdb_group:int) -> list:
        res = []
        df1 = self.heavy_chains(pdb_group)
        g1 = df1.groupby(['pdb_id', 'model_no'])
        df2 = self.light_chains(pdb_group)
        dict2 = self.group_to_dict(df2)
        df3 = self.target_chains(pdb_group)
        dict3 = self.group_to_dict(df3)
        for key, sub1 in g1:
            if key in dict2 and key in dict3:
                sub2, sub3 = dict2[key], dict3[key]
                if len(sub1) == 2 and len(sub2) == 2:
                    chain_combo = pd.concat([sub1['chain_no'], sub2['chain_no']]).to_list()
                    chain_label = pd.concat([sub1['chain_type'], sub2['chain_type']]).to_list()
                    chain_ids = pd.concat([sub1['chain_id'], sub2['chain_id']]).to_list()
                    for _, row3 in sub3.iterrows():
                        res.append({
                            'pdb_id': key[0],
                            'model_no': key[1],
                            'chain_combo': chain_combo + [row3['chain_no'],],
                            'chain_label': chain_label + ['target',],
                            'chain_ids': chain_ids + [row3['chain_id'],],
                        })
        return res
    
    def ab_heavy_seq(self):
        query = f"""
            SELECT DISTINCT chain_id, pdb_id, uniprot_acc, allele_name, pro_seq
            FROM view_antibody_imgt
            WHERE chain_type = 'H'
        ;"""
        rows = self.list_data(query)
        return rows

    def ab_light_seq(self):
        query = f"""
            SELECT chain_id, pdb_id, uniprot_acc, allele_name, pro_seq
            FROM view_antibody_imgt
            WHERE chain_type IN ('L', 'K')
        ;"""
        rows = self.list_data(query)
        return rows

    def ab_seq(self, chain_type:str=None):
        query = f"""
            SELECT chain_id AS id, pdb_id, chain_seq As seq
            FROM view_antibody
            WHERE chain_seq IS NOT NULL
                AND chain_type IN {chain_type};
        """
        rows = self.list_data(query)
        return rows
    
    def ab_uniprot(self):
        query = """
            SELECT * FROM view_pdb_uniprot
            WHERE keyword = 'Immunoglobulin'
                OR keyword LIKE '%antibody%'
        ;"""
        rows = self.list_data(query)
        return rows

    def ab_complex(self):
        query = """
            SELECT * FROM view_pdb_uniprot
            WHERE pdb_id IN(
                SELECT DISTINCT pdb_id FROM view_pdb_uniprot
                WHERE keyword = 'Immunoglobulin'
                    OR keyword LIKE '%antibody%'
            )
        ;"""
        rows = self.list_data(query)
        return rows

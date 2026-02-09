from collections import defaultdict
import os
import pandas as pd

from .query_db import QueryDb

class QueryComplex(QueryDb):

    def __init__(self, verbose:bool=None):
        super().__init__(verbose)
        self.config['database'] = os.getenv('MYSQL_DB', 'COMPLEX2')

    def get_latest_date(self, table_name):
        query = f"""
            SELECT MAX(release_date) AS date FROM {table_name}
        ;"""
        res = self.list_data(query, False)
        return res[0]['date']

    def parse_modelling(self):
        query = """
            SELECT A.chain_id, A.pdb_id, A.uniprot_acc, U.pdb_path
            FROM align_uniprot A, uniprot_structure U
            WHERE A.uniprot_acc = U.uniprot_acc
                AND A.chain_id LIKE "%_0_%"
                AND A.evalue <= 1e-10 AND A.bit_score >= 100
                AND A.chain_id NOT IN (
                    SELECT distinct chain_id
                    FROM uniprot_rmsd
                    WHERE ca_rmsd is NOT NULL
                )
        ;"""
        rows = self.list_data(query)
        return rows
    
    def get_chain_ids(self, pdb_ids:list=None):
        # return all chain ids
        if pdb_ids is None:
            return self.unique_values('chain', 'chain_id')

        chain_ids = []
        # select chain_id belonging to pbd_ids set
        query = """
            SELECT DISTINCT pdb_id, chain_id FROM chain
            WHERE chain_id IS NOT NULL;
        """
        rows = self.list_data(query)
        for row in rows:
            if row['pdb_id'] in pdb_ids:
                if row['chain_id'] not in chain_ids:
                    chain_ids.append(row['chain_id'])
        return chain_ids

    def newly_pdb_ids(self, table_name:str):
        '''
        args: table_name could be compount, chain, keyword, or pdb_method
        pdb_id is foreign key in those tables
        '''
        query = f"""
            SELECT pdb_id, local_path AS pdb_file
            FROM pdb
            WHERE pdb_id NOT IN (
                SELECT DISTINCT pdb_id FROM {table_name}
            );
        """
        rows = self.list_data(query)
        if rows:
             return rows
        return []

    def detect_new_seq(self, df):
        '''
        columns of df is pro_len, seq
        '''
        g = df.groupby('pro_len')
        for pro_len, sub in g:
            input_seq = [i['seq'] for i in sub]
            query = f"""
                SELECT pro_seq AS seq FROM protein_seq
                WHERE pro_len = {pro_len}
            ;"""
            rows = self.list_data(query, True)
            if rows:
                db_seq = [i['seq'] for i in rows]
                new_seq = set(input_seq).difference(db_seq)
                yield pro_len, list(new_seq)

    def parse_pro_id(self, source_name):
        query = f"""
            SELECT pro_id, source_id FROM protein_source
            WHERE source_name = '{source_name}'
        ;"""
        rows = self.list_data(query)
        ids = defaultdict(list)
        for row in rows:
            pro_id = int(row['pro_id'])
            ids[pro_id].append(row['source_id'])
        return ids

    def delete_pdb_ids(self, table_name:str):
        query = f"""
            DELETE FROM {table_name} WHERE chain_type = "x"
        ;"""
        pdb_ids = self.list_data(query)
        return tuple(pdb_ids['pdb_id'])

    def count_complex(self):
        query = f"""
            SELECT complex_type, COUNT(DISTINCT pdb_id) AS pdb_ids
            FROM complex
            GROUP BY complex_type
            ORDER BY complex_type
        ;"""
        res = self.list_data(query)
        return pd.DataFrame(res)
      
    def get_seq(self):
        query = """
            select A.pdb_id, A.compound_no, A.chain_id, A.chain_type, B.pro_seq
            from abag A, chain B
            where A.pdb_id = B.pdb_id
                AND A.compound_no = B.compound_no
                AND A.chain_id = B.chain_id
                AND B.pro_seq is NOT NULL
                AND A.chain_type IN ('light','heavy','single')
        ;"""
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
            INNER JOIN dup     B ON A.pdb_id=B.pdb_id 
                AND A.model_id = B.model_id
                AND A.chain_id = B.chain_id
                AND A.uniprot_acc = B.uniprot_acc
        ;"""
        rows = self.list_data(query)
        print('duplicates of uniprot_align:', len(rows))
        return rows

    def identical_seq_chain(self) -> dict:
        '''
        parse chain_id if their sequences are identical
        '''
        query = "select chain_id, first_chain_id from pdb_chain_seq;"
        df = self.list_data(query, True)

        res = {}
        for _, sub in df.groupby('first_chain_id'):
            ids = sub['chain_id'].tolist()
            for _id in ids:
                res[_id] = ids
        return res
    
    def parse_chain_prediction(self):
        query = """
            SELECT chain_id, predictor, relative_pdb
            FROM chain_plddt
            WHERE relative_pdb REGEXP '_relaxed_'
        """
        pred = self.list_data(query, True)
        g = pred.groupby('chain_id')
        pred = dict(tuple(g))

        res = {}
        equal_chains = self.identical_seq_chain()
        for chain_id, ids in equal_chains.items():
            sub = []
            for equal_id in ids:
                if equal_id in pred:
                    sub += pred[equal_id].to_dict(orient='records')
            if sub:
                res[chain_id] = pd.DataFrame(sub)
        return res
###############################
    def ab_chain_no(self) -> dict:
        query = """
            SELECT pdb_id, chain_no FROM view_antibody
            WHERE model_no = 0
        ;"""
        res = {}
        for row in self.list_data(query):
            pdb_id, chain_no = row['pdb_id'], row['chain_no']
            if pdb_id not in res:
                res[pdb_id] = []
            if chain_no not in res[pdb_id]:
                res[pdb_id].append(chain_no)
        return res

    def ab_combo2(self) -> dict:
        chain_no_pool = self.ab_chain_no()
        res = []
        query = """
            SELECT pdb_id, combo_id, chain_combo
            FROM meta_combo2
            WHERE combo_pdb is not null
        ;"""
        rows = self.list_data(query)
        for row in rows:
            pdb_id = row['pdb_id']
            if pdb_id in chain_no_pool:
                ab_chain_no = set(chain_no_pool[pdb_id])
                combo_chain_no = set(row['chain_combo'].split('-'))
                inter = combo_chain_no.intersection(ab_chain_no)
                if inter:
                    row['ab_chain_no'] = '-'.join(list(inter))
                    res.append(row)
        return pd.DataFrame(res)

    def abag_combo2(self) -> dict:
        chain_no_pool = self.ab_chain_no()
        res = {}
        query = """
            SELECT pdb_id, combo_id, chain_combo
            FROM meta_combo2
            WHERE combo_pdb is not null
        ;"""
        rows = self.list_data(query)
        for row in rows:
            pdb_id = row['pdb_id']
            if pdb_id in chain_no_pool:
                combo_id = row['combo_id']
                ab_chain_no = chain_no_pool[pdb_id]
                chain1, chain2 = row['chain_combo'].split('-')
                # one antibody and one antigen
                if chain1 in ab_chain_no:
                    if chain2 not in ab_chain_no:
                        row['ab_chain_no'] = chain1
                        res[combo_id] = row
                else:
                    if chain2 in ab_chain_no:
                        row['ab_chain_no'] = chain2
                        res[combo_id] = row
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

    def germline_seq(self, chain_type:str):
        query = f"""
            SELECT G.allele_name AS id, S.seq AS seq
            FROM imgt_gene G, imgt_geneseq S
            WHERE G.germline_id = S.germline_id
                AND G.specie = 'Homo_sapiens'
                AND G.isotype = 'IG'
                AND G.chain_type IN {chain_type}
                AND G.region_name = 'V-REGION'
                AND S.seq_type = 'AA-WithoutGaps-F+ORF+inframeP'
        ;"""
        rows = self.list_data(query)
        return rows
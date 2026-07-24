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
        res = self.list_data(query)
        return res[0]['date']

    def get_pdb_groups(self):
        query = "SELECT pdb_group FROM meta_pdb;"
        res = self.list_data(query)
        return [rec['pdb_group'] for rec in res]
    
    def pop_pdb_group(self, table_name):
        query = f"""
            select pdb_group from meta_pdb
            except
            select pdb_group from meta_execute
            where table_name = '{table_name}'
            order by pdb_group
        ;"""
        res = self.list_data(query)
        if len(res) > 0:
            return res[0]['pdb_group']
        return None
    
    def meta_pdb(self, pdb_group:int, table_name:str=None):
        query = f"""
            SELECT pdb_id, relative_pdb
            FROM meta_pdb
            WHERE pdb_group = {pdb_group}
                AND relative_pdb IS NOT NULL"""
        if table_name:
            query += f"""
                AND pdb_id NOT IN (
                    SELECT pdb_id FROM {table_name}
                    WHERE pdb_id IS NOT NULL
                )"""
        rows = self.list_data(query)
        return rows
    
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
            rows = self.list_data(query)
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
        rows = self.list_data(query)
        df = pd.DataFrame(rows)

        res = {}
        for _, sub in df.groupby('first_chain_id'):
            ids = sub['chain_id'].tolist()
            for _id in ids:
                res[_id] = ids
        return res
    
    def parse_chain_prediction(self, table_name:str):
        query = f"""
            SELECT chain_id, relative_pdb
            FROM {table_name}
            WHERE relative_pdb REGEXP '_relaxed_'
        """
        rows = self.list_data(query)
        pred = pd.DataFrame(rows)
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

    def combo_ids(self, pdb_group:int, category:str):
        query = f"""
            SELECT combo_id
            FROM label_complex c
            LEFT JOIN meta_pdb p ON c.pdb_id=p.pdb_id
            WHERE p.pdb_group = {pdb_group}
                AND category= \'{category}\'
        ;"""
        rows = self.list_data(query)
        combo_ids = [i['combo_id'] for i in rows]
        return combo_ids
    
    def pdb_chains(self, pdb_group:int) -> dict:
        query = f"""
            SELECT c.pdb_id, c.model_no, c.chain_id, d.chain_pdb,
                f.relative_path AS chain_faa, s.chain_seq
            FROM pdb_chainid        c
            LEFT JOIN meta_pdb      p ON c.pdb_id=p.pdb_id
            LEFT JOIN pdb_chain_seq s ON c.chain_id=s.chain_id
            LEFT JOIN chain_faa     f ON c.chain_id=f.chain_id
            LEFT JOIN chain_pdb     d ON c.chain_id=d.chain_id
            WHERE p.pdb_group = {pdb_group}
                AND c.model_no = 0
                AND f.relative_path IS NOT NULL
                AND d.chain_pdb IS NOT NULL
        ;"""
        rows = self.list_data(query)
        df = pd.DataFrame(rows)
        g = df.groupby(['pdb_id', 'model_no'])
        res = {key: sub for key, sub in g}
        return res

    def retrieve_combo(self, pdb_group:int, combo:int, dest_table_name:str):
            query = f"""
                SELECT c.pdb_id, c.combo_id, c.combo_pdb
                FROM combo{combo}  c
                LEFT JOIN meta_pdb p ON c.pdb_id = p.pdb_id
                WHERE p.pdb_group = {pdb_group}
                    AND c.combo_pdb IS NOT NULL
                    AND c.combo_id NOT IN (
                        SELECT combo_id FROM {dest_table_name}
                    )
            ;"""
            rows = self.list_data(query)
            return rows

    def pack_proseq(self, query) -> tuple:
        pool = {}
        rows = self.list_data(query)
        for row in rows:
            seq = row['seq']
            seq_head = seq[:10]
            seq_len = len(seq)
            table_name = 'proseq_400' if seq_len <= 400 else ('proseq_1000' \
                if seq_len <=1000 else 'proseq_long')
            row['proseq_name'] = table_name
            key = (table_name, seq_head, seq_len)
            if key not in pool:
                pool[key] = []
            pool[key].append(row)

        #present data
        for key, rows in pool.items():
            if self.verbose:
                print(f'Present {len(key)} data with conditions={key}.')
            yield key, rows

    def parse_seqid(self, conditions, data:list) -> tuple:
        '''
        default key name is 'seq'
        '''
        if not data:
            return ([], [])
        # retrieve seq ~ seq_id
        table_name, seq_head, seq_len = conditions
        qc = QueryComplex(self.verbose)
        query = f"""
            SELECT seq_id, seq FROM {table_name}
            WHERE seq_head = '{seq_head}'
                AND seq_len = {seq_len}
        """
        rows = qc.list_data(query)
        ids = {r['seq']:r['seq_id'] for r in rows}

        #parse, and collect unparsed
        parsed, unknown = [], []
        for item in data:
            seq = item['seq']
            if seq in ids:
                item['seq_id'] = ids[seq]
                parsed.append(item)
            else:
                unknown.append(item)
        return (parsed, unknown)



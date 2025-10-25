import pandas as pd

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

    def get_heavy_seq(self):
        query = f"""
            SELECT DISTINCT chain_id, pdb_id, uniprot_acc, allele_name, pro_seq
            FROM view_antibody_imgt
            WHERE chain_type = 'H';
        """
        rows = self.list_data(query)
        return rows

    def get_light_seq(self):
        query = f"""
            SELECT chain_id, pdb_id, uniprot_acc, allele_name, pro_seq
            FROM view_antibody_imgt
            WHERE chain_type IN ('L', 'K');
        """
        rows = self.list_data(query)
        return rows

    def get_antibody_seq(self, chain_type:str):
        query = f"""
            SELECT chain_id AS id, pdb_id, pro_seq As seq
            FROM view_antibody_imgt
            WHERE chain_type IN {chain_type};
        """
        rows = self.list_data(query)
        return rows

    def get_germline_seq(self, chain_type:str):
        query = f"""
            SELECT G.allele_name AS id, S.seq AS seq
            FROM imgt_gene G, imgt_geneseq S
            WHERE G.germline_id = S.germline_id
                AND G.specie = 'Homo_sapiens'
                AND G.isotype = 'IG'
                AND G.chain_type IN {chain_type}
                AND G.region_name = 'V-REGION'
                AND S.seq_type = 'AA-WithoutGaps-F+ORF+inframeP';
        """
        rows = self.list_data(query)
        return rows

    def germline_human_vregion(self, pdb_ids:list):
        records_pool = {}
        pool = [
            ["H", "('H')"],
            ["K", "('K')"],
            ["L", "('L')"],
            ["KL", "('K', 'L')"],
        ]
        for chain_type, cond in pool:
            records_pool[chain_type] = []
            germlines = self.get_germline_seq(cond)
            abs = self.get_antibody_seq(cond)
            for ab in abs:
                if ab['pdb_id'] in pdb_ids:
                    item = (ab, germlines)
                    records_pool[chain_type].append(item)
        return records_pool
    
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
                );
        """
        rows = self.list_data(query)
        return rows
    

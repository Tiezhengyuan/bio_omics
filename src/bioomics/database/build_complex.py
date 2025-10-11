'''
insert data into empty table
'''
import os
import numpy as np
import pandas as pd

# from src.process_pdb import ProcessPdb
from .build_db import BuildDb
from bioomics import Pdb, PullPdb, ParseDatRecord, ParseBlastXml

class BuildComplex(BuildDb):

    def __init__(self, data:list=None, verbose:bool=False):
        '''
        args: data is iterable of pdb_file, or dat_record
        '''
        super().__init__(data, verbose)
        self.config = {
            'host': 'localhost',
            'user': 'admin',
            'password': 'strong_password',
            'database': 'COMPLEX2'
        }
    
    def record_pdb(self, pdb_file):
        p = Pdb(pdb_file, None, self.verbose)
        p.load_structure()
        pdb_id = p.structure_id

        pdb_url = f"https://www.rcsb.org/structure/{pdb_id}"
        release_date = p.structure.header['release_date']
        name = p.structure.header['name']
        resolution = p.structure.header['resolution']
        rec = (pdb_id, pdb_url, pdb_file, release_date, resolution, name)
        return [rec,]

    def insert_pdb(self):
        # insert data into table pdb
        table_name = 'pdb'
        insert_query = """
            INSERT INTO pdb (pdb_id, pdb_url, local_path, release_date, resolution, name)
            VALUES (%s, %s, %s, %s, %s, %s)
        """
        self.insert(table_name, insert_query, self.record_pdb)

    def record_compound(self, pdb_file):
        # records
        p = Pdb(pdb_file, None, self.verbose)
        p.load_structure()
        pdb_id = p.structure_id

        compound = p.structure.header['compound']
        records = []
        for compound_no, cp in compound.items():
            compound_no = int(compound_no)
            chain = cp.get('chain').replace(' ', '')
            _text = []
            for _key, _value in cp.items():
                if _key not in ('chain',) and _value:
                    _text.append(f"{_key}: {_value}")
            compound_text = ' ; '.join(_text)
            rec = (pdb_id, compound_no, chain, compound_text)
            records.append(rec)
        return records

    def insert_compound(self) -> tuple:
        table_name = 'compound'
        insert_query = """
            INSERT INTO compound (pdb_id, compound_no, chain, compound_text)
            VALUES (%s, %s, %s, %s)
        """
        self.insert(table_name, insert_query, self.record_compound)

    def record_chain(self, pdb_file):
        p = Pdb(pdb_file, None, self.verbose)
        p.load_structure()
        chains = p.get_chains_seq()
        pdb_id = p.structure_id

        records = []
        compound = p.structure.header['compound']
        for compound_no, cp in compound.items():
            compound_no = int(compound_no)
            cp_chains = cp.get('chain').replace(' ', '').split(',')
            for chain_id in cp_chains:
                chain_id = chain_id.upper()
                pro_seq = chains.get(chain_id)
                rec = (pdb_id, compound_no, chain_id, pro_seq)
                if rec not in records:
                    records.append(rec)
        return records

    def insert_chain(self) -> tuple:
        table_name = 'chain'
        insert_query = """
            INSERT INTO chain (pdb_id, compound_no, chain_id, pro_seq)
            VALUES (%s, %s, %s, %s)
        """
        self.insert(table_name, insert_query, self.record_chain)

    def record_pdb_keyword(self, pdb_file):
        p = Pdb(pdb_file, None, self.verbose)
        p.load_structure()
        pdb_id = p.structure_id

        records = []
        keywords = p.structure.header['keywords']
        keywords = keywords.split(', ')
        for keyword in keywords:
            rec = (pdb_id, keyword)
            records.append(rec)
        return records

    def insert_pdb_keyword(self) -> tuple:
        table_name = 'pdb_keyword'
        insert_query = """
            INSERT INTO pdb_keyword (pdb_id, keyword)
            VALUES (%s, %s)
        """
        self.insert(table_name, insert_query, self.record_pdb_keyword, 500)

    def record_imgt_genelist(self, record:dict):
        names = ['Species', 'IMGT/GENE-DB', 'IMGT gene functionality', \
            'IMGT and HGNC gene definition', 'Number of alleles', \
            'Chromosome', 'Chromosomal localization', \
            'IMGT/LIGM-DB reference sequence(s) for allele *01', \
            'NCBI Gene', 'chain', 'domain']
        rec = [record.get(i, '') for i in names]
        return [rec,]

    def insert_imgt_genelist(self) -> tuple:
        table_name = 'imgt_genelist'
        insert_query = f"""
            INSERT INTO {table_name} (species, gene_db, gene_function, \
                gene_definition, number_alleles, chromosome, chromosomal_loc, \
                ligm_db_ref_gene, ncbi_gene, chain, domain)
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
        """
        self.insert(table_name, insert_query, self.record_imgt_genelist)

########################
    def record_uniprot(self, record):
        p = ParseDatRecord(record)
        rec = (p.uniprot_acc, p.name, p.organism, p.seq,)
        return [rec,]

    def insert_uniprot(self) -> tuple:
        table_name = 'uniprot'
        insert_query = f"""
            INSERT INTO {table_name} (uniprot_acc, name, organism, pro_seq)
            VALUES (%s, %s, %s, %s)
        """
        self.insert(table_name, insert_query, self.record_uniprot)

    def record_uniprot_keyword(self, record):
        uniprot_acc = record.id
        records = []
        keywords = ParseDatRecord(record).keywords()
        for keyword in keywords:
            rec = (uniprot_acc, keyword)
            if rec not in records:
                records.append(rec)
        return records

    def insert_uniprot_keyword(self) -> tuple:
        table_name = 'uniprot_keyword'
        insert_query = f"""
            INSERT INTO {table_name} (uniprot_acc, keyword)
            VALUES (%s, %s)
        """
        self.insert(table_name, insert_query, self.record_uniprot_keyword, 500)

    def record_uniprot_dbxrefs(self, record):
        uniprot_acc = record.id
        records = []
        dbxrefs = ParseDatRecord(record).dbxrefs()
        for k, values in dbxrefs.items():
            for v in values:
                rec = (uniprot_acc, k, v)
                if rec not in records:
                    records.append(rec)
        return records

    def insert_uniprot_dbxrefs(self) -> tuple:
        table_name = 'uniprot_dbxrefs'
        insert_query = f"""
            INSERT INTO {table_name} (uniprot_acc, db_name, db_acc)
            VALUES (%s, %s, %s)
        """
        self.insert(table_name, insert_query, self.record_uniprot_dbxrefs, 1000)


    def record_uniprot_align(self, aln_file):
        best_hits = ParseBlastXml(aln_file).best_hits()
        records = []
        for hit in best_hits:
            pdb_id, model_id, chain_id = hit['query'].split('_')
            uniprot_acc = hit['hit_desc'].split('|')[1]
            rec = (pdb_id, int(model_id), chain_id.upper(), \
                uniprot_acc, float(hit['evalue']), \
                int(hit['bit_score']), float(hit['identity']))
            if rec not in records:
                records.append(rec)
        return records

    def insert_uniprot_align(self) -> tuple:
        table_name = 'uniprot_align'
        insert_query = f"""
            INSERT INTO {table_name} (pdb_id, model_id, chain_id, 
                uniprot_acc, evalue, bit_score, identity)
            VALUES (%s, %s, %s, %s, %s, %s, %s)
        """
        self.insert(table_name, insert_query, self.record_uniprot_align)


#####################

    def insert_complex(self):
        '''
        insert data into empty table complex
        '''
        cols = ['pdb_id', 'complex_type']
        df = self.data[cols].drop_duplicates()
        data = df.to_records(index=False)
        insert_query = """
            INSERT INTO complex (pdb_id, complex_type)
            VALUES (%s, %s)
        """
        self.insert_data(insert_query, data)

    def insert_abag(self):
        '''
        table abag
        '''
        # avoid error of int64
        self.data['compound_no'] = self.data['compound_no'].astype(object)
        df = self.data.where(pd.notnull(self.data), None)
        # column: chain_id
        df = df.assign(chain_id=df['chain'].str.split(',')).explode('chain_id')
        df['chain_id'] = df['chain_id'].str.upper()

        # select
        cols = ['pdb_id', 'compound_no', 'chain_id', 'chain_type', 'fragment']
        data = df[cols].to_records(index=False)
        insert_query = """
            INSERT INTO abag(pdb_id, compound_no, chain_id, chain_type, fragment)
            VALUES (%s, %s, %s, %s, %s)
        """
        self.insert_data(insert_query, data)


    def record_aacdb(self):
        '''
        self.data is df of aacdb summary
        '''
        for ix, row in self.data.iterrows():
            pdb_id = row['pdb'].upper()
            chains = row['chains'].split('_')
            rec = (
                pdb_id,
                chains[0],
                chains[1] if len(chains) > 1 else None,
                row['antibody'] ,
                row['protein'],
                None if np.isnan(row['resolution']) else row['resolution'],
            )
            yield pdb_id, [rec,]


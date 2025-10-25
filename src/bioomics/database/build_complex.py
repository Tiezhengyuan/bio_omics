'''
insert data into empty table
'''
import numpy as np
import pandas as pd


from .build_db import BuildDb
from .query_complex import QueryComplex
from ..PDB import Pdb
from ..UniProt import ParseDatRecord
from ..protein import ParseBlastXml, ParseIgblast

class BuildComplex(BuildDb):

    def __init__(self, data:list=None, verbose:bool=False, chunk_size:int=None):
        '''
        args: data is iterable of pdb_file, or dat_record
        '''
        super().__init__(data, verbose, chunk_size)
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
        insert_query = f"""INSERT INTO {table_name}
            (pdb_id, pdb_url, local_path, release_date, resolution, name)
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
            chains = cp.get('chain').replace(' ', '').upper()
            _text = []
            for _key, _value in cp.items():
                if _key not in ('chain',) and _value:
                    _text.append(f"{_key}: {_value}")
            compound_text = ' ; '.join(_text)
            rec = (pdb_id, compound_no, chains, compound_text)
            records.append(rec)
        return records

    def insert_compound(self) -> tuple:
        table_name = 'compound'
        insert_query = f"""INSERT INTO {table_name} 
            (pdb_id, compound_no, chains, compound_text)
            VALUES (%s, %s, %s, %s)"""
        self.insert(table_name, insert_query, self.record_compound)

    def record_chain(self, pdb_file):
        p = Pdb(pdb_file, None, self.verbose)
        p.load_structure()
        pdb_id = p.structure_id
        records = []
        for model_chains in p.get_chains():
            model_no = int(model_chains['model_no'])
            for chain in model_chains['chains']:
                chain_no = chain['chain_no']
                chain_id = f"{pdb_id}_{model_no}_{chain_no}"
                pro_seq = chain['seq']
                # the column pro_seq is not null
                if pro_seq:
                    rec = (chain_id, pdb_id, model_no, chain_no, pro_seq)
                    if rec not in records:
                        records.append(rec)
        return records

    def record_check_chain(self, pdb_file:str, chain_ids:list):
        checked_records = []
        for rec in self.record_chain(pdb_file):
            chain_id = rec[0]
            if chain_id not in chain_ids:
                checked_records.append(rec)
        return checked_records

    def insert_chain(self, check:bool=False) -> tuple:
        table_name = 'chain'
        insert_query = f"""INSERT INTO {table_name} 
            (chain_id, pdb_id, model_no, chain_no, pro_seq)
            VALUES (%s, %s, %s, %s, %s)"""
        if check:
            chain_ids = QueryComplex().get_unique_values(table_name, 'chain_id')
            return self.insert(table_name, insert_query, self.record_check_chain, \
                chain_ids=chain_ids)
        return self.insert(table_name, insert_query, self.record_chain)

    def record_keyword(self, pdb_file):
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

    def insert_keyword(self) -> tuple:
        table_name = 'keyword'
        insert_query = f"""INSERT INTO {table_name} 
            (pdb_id, keyword) VALUES (%s, %s)"""
        self.insert(table_name, insert_query, self.record_keyword)

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
        insert_query = f"""INSERT INTO {table_name} (species, gene_db, 
                gene_function, gene_definition, number_alleles, chromosome,
                chromosomal_loc, ligm_db_ref_gene, ncbi_gene, chain, domain)
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)"""
        self.insert(table_name, insert_query, self.record_imgt_genelist)

    def record_imgt_gene(self, record:dict):
        names = ['germline_id', 'specie', 'isotype', 'chain_type', \
            'gene_name', 'allele_name', 'region_name']
        rec = [record.get(i, '') for i in names]
        return [rec,]

    def insert_imgt_gene(self) -> tuple:
        table_name = 'imgt_gene'
        insert_query = f"""INSERT INTO {table_name} (germline_id, specie,
                isotype, chain_type, gene_name, allele_name, region_name)
            VALUES (%s, %s, %s, %s, %s, %s, %s)"""
        self.insert(table_name, insert_query, self.record_imgt_gene)

    def record_imgt_geneseq(self, record:dict):
        names = ['germline_id', 'seq_type', 'seq']
        rec = [record.get(i, '') for i in names]
        return [rec,]

    def insert_imgt_geneseq(self) -> tuple:
        table_name = 'imgt_geneseq'
        insert_query = f"""INSERT INTO {table_name}
            (germline_id, seq_type, seq) VALUES (%s, %s, %s)"""
        self.insert(table_name, insert_query, self.record_imgt_geneseq)

    def record_uniprot_structure(self, record:dict):
        names = ['uniprot_acc', 'model', 'pdb_path']
        rec = [record.get(i, '') for i in names]
        return [rec,]

    def insert_uniprot_structure(self) -> tuple:
        table_name = 'uniprot_structure'
        insert_query = f"""INSERT INTO {table_name}
            (uniprot_acc, model, pdb_path) VALUES (%s, %s, %s)"""
        self.insert(table_name, insert_query, self.record_uniprot_structure)

########################
    def record_uniprot(self, record):
        p = ParseDatRecord(record)
        rec = (p.uniprot_acc, p.name, p.organism, p.seq,)
        return [rec,]

    def insert_uniprot(self) -> tuple:
        table_name = 'uniprot'
        insert_query = f"""INSERT INTO {table_name}
            (uniprot_acc, name, organism, pro_seq)
            VALUES (%s, %s, %s, %s)"""
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
        insert_query = f"""INSERT INTO {table_name}
            (uniprot_acc, keyword) VALUES (%s, %s)"""
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
        insert_query = f"""INSERT INTO {table_name}
            (uniprot_acc, db_name, db_acc)
            VALUES (%s, %s, %s)"""
        self.insert(table_name, insert_query, self.record_uniprot_dbxrefs)

    def record_align_uniprot(self, aln_file):
        records = []
        best_hits = ParseBlastXml(aln_file).best_hits()
        for hit in best_hits:
            chain_id = hit['query'].upper()
            pdb_id = hit['query'].split('_')[0]
            rec0 = (chain_id, pdb_id)
            if 'hit_desc' in hit:
                uniprot_acc = hit['hit_desc'].split('|')[1]
                rec = rec0 + (
                    uniprot_acc,
                    float(hit['evalue']),
                    int(hit['bit_score']),
                    float(hit['identity'])
                )
                if rec not in records:
                    records.append(rec)
            else:
                # insert null if not hits detected
                rec0 += (None,) * 4
                if rec0 not in records:
                    records.append(rec0)
        return records

    def insert_align_uniprot(self) -> tuple:
        table_name = 'align_uniprot'
        insert_query = f"""INSERT INTO {table_name}
            (chain_id, pdb_id, uniprot_acc, evalue, bit_score, identity)
            VALUES (%s, %s, %s, %s, %s, %s)"""
        self.insert(table_name, insert_query, self.record_align_uniprot)

    def record_align_imgt(self, item:tuple):
        aln_file, fa_file = item
        pi = ParseIgblast(aln_file, fa_file, self.verbose)
        records = []
        for hit in pi.retrieve_hits():
            chain_id = hit['query_id']
            pdb_id = chain_id.split('_')[0]
            rec0 = (chain_id, pdb_id)
            if 'allele_name' in hit:
                rec = rec0 + (
                    hit['allele_name'].split('|')[0],
                    hit['gene_name'],
                    hit['gene_family'],
                    hit['chain_type'].upper(),
                    hit['isotype'].upper(), 
                    float(hit['e-value']),
                    float(hit['bits-score']),
                )
                if rec not in records:
                    records.append(rec)
            else:
                # insert null if not hits detected
                rec0 += (None,) * 7
                if rec0 not in records:
                    records.append(rec0)
        return records

    def insert_align_imgt(self) -> tuple:
        table_name = 'align_imgt'
        insert_query = f"""INSERT INTO {table_name} (chain_id,
            pdb_id, allele_name, gene_name, gene_family,
            chain_type, isotype, evalue, bit_score)
            VALUES (%s, %s, %s, %s, %s, %s,%s, %s, %s)
        """
        self.insert(table_name, insert_query, self.record_align_imgt)

    def record_align_vregion(self, item:tuple):
        aln_file, fa_file = item
        pi = ParseIgblast(aln_file, fa_file, self.verbose)
        records = []
        for region in pi.retrieve_regions():
            chain_id = region['query_id']
            pdb_id = chain_id.split('_')[0]
            rec0 = (chain_id, pdb_id,)
            if 'percent identity' in region:
                rec = rec0 + (
                    region['region'],
                    region['percent identity'],
                    region['from'],
                    region['to'],
                    region['seq'],
                )
                if rec not in records:
                    records.append(rec)
            else:
                # region name is composite primary key
                rec0 += ('UN', None, None, None, None)
                if rec0 not in records:
                    records.append(rec0)
        return records

    def insert_align_vregion(self) -> tuple:
        table_name = 'align_vregion'
        insert_query = f"""INSERT INTO {table_name}
            (chain_id, pdb_id, region_name, identity, start, end, seq)
            VALUES (%s, %s, %s, %s, %s, %s, %s)
        """
        self.insert(table_name, insert_query, self.record_align_vregion)

    def insert_uniprot_rmsd(self, record):
        table_name = 'uniprot_rmsd'
        insert_query = f"""INSERT INTO {table_name}
            (chain_id, native_pdb, uniprot_acc, model_pdb, ca_rmsd)
            VALUES (%s, %s, %s, %s, %s)
        """
        records = [record,]
        self.insert_records(insert_query, records)

#####################

    def insert_complex(self):
        '''
        insert data into empty table complex
        '''
        cols = ['pdb_id', 'complex_type']
        df = self.data[cols].drop_duplicates()
        data = df.to_records(index=False)
        insert_query = """INSERT INTO complex
            (pdb_id, complex_type) VALUES (%s, %s)"""
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
        insert_query = """INSERT INTO abag
            (pdb_id, compound_no, chain_id, chain_type, fragment)
            VALUES (%s, %s, %s, %s, %s)"""
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


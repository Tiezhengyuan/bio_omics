'''
insert data into empty table
'''
import numpy as np
import pandas as pd


from .build_db import BuildDb
from .query_complex import QueryComplex
from .build_records import BuildRecords

from ..protein import ParseIgblast


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
        self.recorder = BuildRecords(verbose)
    
    def insert_pdb(self):
        # insert data into table pdb
        table_name = 'pdb'
        insert_query = f"""INSERT INTO {table_name}
            (pdb_id, pdb_url, release_date, name)
            VALUES (%s, %s, %s, %s)"""
        func = getattr(self.recorder, table_name)
        self.insert(table_name, insert_query, func)


    def insert_pdb_method(self):
        # insert data into table pdb
        table_name = 'pdb_method'
        insert_query = f"""INSERT INTO {table_name}
            (pdb_id, structure_method, resolution)
            VALUES (%s, %s, %s)"""
        func = getattr(self.recorder, table_name)
        self.insert(table_name, insert_query, func)

    def insert_pdb_compound(self) -> tuple:
        table_name = 'pdb_compound'
        insert_query = f"""INSERT INTO {table_name} 
            (pdb_id, compound_no, chains, compound_text)
            VALUES (%s, %s, %s, %s)"""
        func = getattr(self.recorder, table_name)
        self.insert(table_name, insert_query, func)

    def insert_pdb_keyword(self) -> tuple:
        table_name = 'pdb_keyword'
        insert_query = f"""INSERT INTO {table_name} 
            (pdb_id, keyword) VALUES (%s, %s)"""
        func = getattr(self.recorder, table_name)
        self.insert(table_name, insert_query, func)

    def insert_pdb_chainid(self) -> tuple:
        table_name = 'pdb_chainid'
        insert_query = f"""INSERT INTO {table_name} 
            (pdb_id, model_no, chain_no, chain_id)
            VALUES (%s, %s, %s, %s)"""
        func = getattr(self.recorder, table_name)
        return self.insert(table_name, insert_query, func)
    
    def insert_pdb_chain(self) -> tuple:
        table_name = 'pdb_chain'
        insert_query = f"""INSERT INTO {table_name} 
            (chain_id, pdb_id, model_no, chain_no, chain_seq)
            VALUES (%s, %s, %s, %s, %s)"""
        func = getattr(self.recorder, table_name)
        return self.insert(table_name, insert_query, func)

    def insert_pdb_chain_seq(self) -> tuple:
        table_name = 'pdb_chain_seq'
        insert_query = f"""INSERT INTO {table_name} 
            (chain_id, chain_seq) VALUES (%s, %s)"""
        func = getattr(self.recorder, table_name)
        return self.insert(table_name, insert_query, func)

    def insert_meta_bfactor(self) -> tuple:
        table_name = 'meta_bfactor'
        insert_query = f"""INSERT INTO {table_name} 
            (pdb_id, pickle_file) VALUES (%s, %s)"""
        func = lambda rec: [rec,]
        self.insert(table_name, insert_query, func)

    def insert_pdb_bfactor(self) -> tuple:
        table_name = 'pdb_bfactor'
        insert_query = f"""INSERT INTO {table_name} 
            (pdb_id, model_no, chain_no, avg_bfactor)
            VALUES (%s, %s, %s, %s)"""
        func = getattr(self.recorder, table_name)
        self.insert(table_name, insert_query, func)

    def insert_pdb_secondary(self) -> tuple:
        table_name = 'pdb_secondary'
        insert_query = f"""INSERT INTO {table_name} 
            (pdb_id, model_no, chain_no, dssp_seq)
            VALUES (%s, %s, %s, %s)"""
        func = getattr(self.recorder, table_name)
        self.insert(table_name, insert_query, func)

    def insert_align_uniprot(self) -> tuple:
        table_name = 'align_uniprot'
        insert_query = f"""INSERT INTO {table_name}
            (chain_id, pdb_id, uniprot_acc, evalue, bit_score, identity)
            VALUES (%s, %s, %s, %s, %s, %s)"""
        func = getattr(self.recorder, table_name)
        self.insert(table_name, insert_query, func, 2)

########################
    def insert_meta_pdb(self):
        table_name = 'meta_pdb'
        insert_query = f"""INSERT INTO {table_name}
            (pdb_id, pdb_file, pdb_group)
            VALUES (%s, %s, %s)"""
        func = lambda rec: [rec,]
        self.insert(table_name, insert_query, func)

    def insert_meta_faa(self):
        table_name = 'meta_faa'
        insert_query = f"""INSERT INTO {table_name}
            (pdb_id, faa_file) VALUES (%s, %s)"""
        func = lambda rec: [rec,]
        self.insert(table_name, insert_query, func)

    def insert_meta_chain_faa(self):
        table_name = 'meta_chain_faa'
        insert_query = f"""INSERT INTO {table_name}
            (chain_id, faa_file) VALUES (%s, %s)"""
        func = lambda rec: [rec,]
        self.insert(table_name, insert_query, func)

    def insert_meta_chain_pdb(self):
        table_name = 'meta_chain_pdb'
        # chain_no is case sensitive in PDB
        insert_query = f"""INSERT INTO {table_name}
            (chain_id, pdb_file) VALUES (%s, %s)"""
        func = lambda rec: [rec,]
        self.insert(table_name, insert_query, func)

    def insert_meta_combo_pdb(self):
        table_name = 'meta_combo_pdb'
        # chain_no is case sensitive in PDB
        insert_query = f"""INSERT INTO {table_name}
            (pdb_id, model_no, chain_combo, pdb_file)
            VALUES (%s, %s, %s, %s)"""
        func = lambda rec: [rec,]
        self.insert(table_name, insert_query, func)

    def insert_meta_blastp5(self):
        table_name = 'meta_blastp5'
        # chain_no is case sensitive in PDB
        insert_query = f"""INSERT INTO {table_name}
            (pdb_id, aln_file) VALUES (%s, %s)"""
        func = lambda rec: [rec,]
        self.insert(table_name, insert_query, func)

    def insert_meta_pdb_igblastp(self):
        table_name = 'meta_pdb_igblastp'
        # chain_no is case sensitive in PDB
        insert_query = f"""INSERT INTO {table_name}
            (pdb_id, aln_file) VALUES (%s, %s)"""
        func = lambda rec: [rec,]
        self.insert(table_name, insert_query, func)

    def insert_meta_pdb_distance(self):
        table_name = 'meta_pdb_distance'
        # chain_no is case sensitive in PDB
        insert_query = f"""INSERT INTO {table_name}
            (pdb_id, pickle_file) VALUES (%s, %s)"""
        func = lambda rec: [rec,]
        self.insert(table_name, insert_query, func)

    def insert_meta_pdb_sasa(self):
        table_name = 'meta_pdb_sasa'
        # chain_no is case sensitive in PDB
        insert_query = f"""INSERT INTO {table_name}
            (pdb_id, pickle_file) VALUES (%s, %s)"""
        func = lambda rec: [rec,]
        self.insert(table_name, insert_query, func)

######################
    def insert_imgt_genelist(self) -> tuple:
        table_name = 'imgt_genelist'
        insert_query = f"""INSERT INTO {table_name}
            (species, gene_db, gene_function, gene_definition,
                number_alleles, chromosome, chromosomal_loc,
                ligm_db_ref_gene, ncbi_gene, chain, domain)
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)"""
        func = getattr(self.recorder, table_name)
        self.insert(table_name, insert_query, func)

    def insert_imgt_gene(self) -> tuple:
        table_name = 'imgt_gene'
        insert_query = f"""INSERT INTO {table_name}
            (germline_id, specie, isotype, chain_type,
                gene_name, allele_name, region_name)
            VALUES (%s, %s, %s, %s, %s, %s, %s)"""
        func = getattr(self.recorder, table_name)
        self.insert(table_name, insert_query, func)

    def insert_imgt_geneseq(self) -> tuple:
        table_name = 'imgt_geneseq'
        insert_query = f"""INSERT INTO {table_name}
            (germline_id, seq_type, seq)
            VALUES (%s, %s, %s)"""
        func = getattr(self.recorder, table_name)    
        self.insert(table_name, insert_query, func)

    def insert_msa_vregion(self) -> tuple:
        table_name = 'msa_vregion'
        insert_query = f"""INSERT INTO {table_name}
            (chain_id, chain_types, query_fa, aln_file, pickle_file)
            VALUES (%s, %s, %s, %s, %s)"""
        func = lambda rec: [rec,]
        self.insert(table_name, insert_query, func)

########################

    def insert_uniprot(self) -> tuple:
        table_name = 'uniprot'
        insert_query = f"""INSERT INTO {table_name}
            (uniprot_acc, name, organism, pro_seq)
            VALUES (%s, %s, %s, %s)"""
        func = getattr(self.recorder, table_name)
        self.insert(table_name, insert_query, func)

    def insert_uniprot_keyword(self) -> tuple:
        table_name = 'uniprot_keyword'
        insert_query = f"""INSERT INTO {table_name}
            (uniprot_acc, keyword) VALUES (%s, %s)"""
        func = getattr(self.recorder, table_name)
        self.insert(table_name, insert_query, func)

    def insert_uniprot_dbxrefs(self) -> tuple:
        table_name = 'uniprot_dbxrefs'
        insert_query = f"""INSERT INTO {table_name}
            (uniprot_acc, db_name, db_acc)
            VALUES (%s, %s, %s)"""
        func = getattr(self.recorder, table_name)
        self.insert(table_name, insert_query, func)

    def insert_uniprot_meta_alphafold(self) -> tuple:
        table_name = 'uniprot_meta_alphafold'
        insert_query = f"""INSERT INTO {table_name}
            (uniprot_acc, pdb_path) VALUES (%s, %s)"""
        func = lambda rec: [rec,]
        self.insert(table_name, insert_query, func)

    def insert_uniprot_alphafold_plddt(self) -> tuple:
        table_name = 'uniprot_alphafold_plddt'
        insert_query = f"""INSERT INTO {table_name}
            (uniprot_acc, model_no, chain_no, avg_plddt, pickle_file)
            VALUES (%s, %s, %s, %s, %s)"""
        func = getattr(self.recorder, table_name)
        self.insert(table_name, insert_query, func)

    def insert_align_imgt(self) -> tuple:
        table_name = 'align_imgt'
        insert_query = f"""INSERT INTO {table_name} (chain_id,
            pdb_id, allele_name, gene_name, gene_family,
            chain_type, isotype, evalue, bit_score)
            VALUES (%s, %s, %s, %s, %s, %s,%s, %s, %s)
        """
        func = getattr(self.recorder, table_name)
        self.insert(table_name, insert_query, func)

    def insert_align_vregion(self) -> tuple:
        table_name = 'align_vregion'
        insert_query = f"""INSERT INTO {table_name}
            (chain_id, pdb_id, region_name, identity, start, end, seq)
            VALUES (%s, %s, %s, %s, %s, %s, %s)
        """
        func = getattr(self.recorder, table_name)
        self.insert(table_name, insert_query, func)

    def insert_uniprot_alphafold_rmsd(self):
        table_name = 'uniprot_alphafold_rmsd'
        insert_query = f"""INSERT INTO {table_name}
            (chain_id, uniprot_acc, ca_rmsd, native_pdb, model_pdb)
            VALUES (%s, %s, %s, %s, %s)
        """
        func = lambda rec: [rec,]
        self.insert(table_name, insert_query, func)

    def insert_uniprot_rmsd(self):
        table_name = 'uniprot_rmsd'
        insert_query = f"""INSERT INTO {table_name}
            (chain_id, native_pdb, uniprot_acc, model_pdb, ca_rmsd)
            VALUES (%s, %s, %s, %s, %s)
        """
        func = lambda rec: [rec,]
        self.insert(table_name, insert_query, func)

#####################


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

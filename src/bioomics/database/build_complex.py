'''
insert data into empty table
'''
import os
from .build_db import BuildDb
from .build_records import BuildRecords


class BuildComplex(BuildDb):

    def __init__(self, data:list=None, verbose:bool=False, chunk_size:int=None):
        '''
        args: data is iterable of pdb_file, or dat_record
        '''
        super().__init__(data, verbose, chunk_size)
        self.config['database'] = os.getenv('MYSQL_DB', 'COMPLEX2')
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
            (chain_id, pro_len, chain_seq)
            VALUES (%s, %s, %s)"""
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

    def insert_pdb_chain_type(self) -> tuple:
        table_name = 'pdb_chain_type'
        insert_query = f"""INSERT INTO {table_name} 
            (pdb_id, chain_id, chain_type, total_res, non_water_res,
                aa_res, dna_res, rna_res, water_res, ion_res,
                hetero_res, other_res)
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
        """
        func = lambda rec: [rec,]
        return self.insert(table_name, insert_query, func)
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

    def insert_meta_pdb_fixed(self):
        table_name = 'meta_pdb_fixed'
        insert_query = f"""INSERT INTO {table_name}
            (pdb_id, pdb_file, pdb_group)
            VALUES (%s, %s, %s)"""
        func = lambda rec: [rec,]
        self.insert(table_name, insert_query, func)

    def insert_meta_pdb_relaxed(self):
        table_name = 'meta_pdb_relaxed'
        insert_query = f"""INSERT INTO {table_name}
            (pdb_id, pdb_file, pdb_group)
            VALUES (%s, %s, %s)"""
        func = lambda rec: [rec,]
        self.insert(table_name, insert_query, func)

######################
    def insert_imgt_genelist(self) -> tuple:
        table_name = 'imgt_genelist'
        insert_query = f"""INSERT INTO {table_name}
            (species, gene_db, gene_function, gene_definition,
                number_alleles, chromosome, chromosomal_loc,
                ligm_db_ref_gene, ncbi_gene, chain, domain)
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
        """
        func = getattr(self.recorder, table_name)
        self.insert(table_name, insert_query, func)

    def insert_imgt_gene(self) -> tuple:
        table_name = 'imgt_gene'
        insert_query = f"""INSERT INTO {table_name}
            (record_id, specie, isotype, chain_type,
                gene_name, allele_name, region_name)
            VALUES (%s, %s, %s, %s, %s, %s, %s)"""
        func = lambda rec: [rec,]
        self.insert(table_name, insert_query, func)

    def insert_imgt_geneseq(self) -> tuple:
        table_name = 'imgt_geneseq'
        insert_query = f"""INSERT INTO {table_name}
            (record_id, seq_type, seq) VALUES (%s, %s, %s)"""
        func = lambda rec: [rec,]
        self.insert(table_name, insert_query, func)

    def insert_msa_vregion(self) -> tuple:
        table_name = 'msa_vregion'
        insert_query = f"""INSERT INTO {table_name}
            (chain_id, chain_types, query_fa, aln_file, pickle_file)
            VALUES (%s, %s, %s, %s, %s)"""
        func = lambda rec: [rec,]
        self.insert(table_name, insert_query, func)

    def insert_imgt_inn(self) -> tuple:
        table_name = 'imgt_inn'
        insert_query = f"""INSERT INTO {table_name}
            (inn_file, inn_number, cas_number, common_name, 
                inn_name, proposed_list, recommended_list, 
                receptor_type, receptor_description, species,
                molecular_formula, glycosylation_sites)
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
        """
        func = getattr(self.recorder, table_name)
        self.insert(table_name, insert_query, func)

    def insert_imgt_inn_chain(self) -> tuple:
        table_name = 'imgt_inn_chain'
        insert_query = f"""INSERT INTO {table_name}
            (inn_number, chain_id, chain_description, pro_len, chain_seq)
            VALUES (%s, %s, %s, %s, %s)
        """
        func = getattr(self.recorder, table_name)
        self.insert(table_name, insert_query, func)

    def insert_imgt_inn_region(self) -> tuple:
        table_name = 'imgt_inn_region'
        insert_query = f"""INSERT INTO {table_name}
            (chain_id, region_name, region_type, start, end, seq)
            VALUES (%s, %s, %s, %s, %s, %s)
        """
        func = getattr(self.recorder, table_name)
        self.insert(table_name, insert_query, func)

    def insert_imgt_inn_cdomain(self) -> tuple:
        table_name = 'imgt_inn_cdomain'
        insert_query = f"""INSERT INTO {table_name}
            (chain_id, domain_name, gene_allele, seq, seq_align)
            VALUES (%s, %s, %s, %s, %s)
        """
        func = getattr(self.recorder, table_name)
        self.insert(table_name, insert_query, func)

    def insert_imgt_inn_vdomain(self) -> tuple:
        table_name = 'imgt_inn_vdomain'
        insert_query = f"""INSERT INTO {table_name}
            (chain_id, domain_name, cdr_length, fr_length,
            gene_allele, seq, seq_cdr, seq_align, seq_align_cdr,
            seq_cdr_mask, seq_align_cdr_mask)
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
        """
        func = getattr(self.recorder, table_name)
        self.insert(table_name, insert_query, func)

    def insert_imgt_inn_cdr(self) -> tuple:
        table_name = 'imgt_inn_cdr'
        insert_query = f"""INSERT INTO {table_name}
            (chain_id, cdr_type, seq, seq_align,
            align_start, align_end, align_length)
            VALUES (%s, %s, %s, %s, %s, %s, %s)
        """
        func = getattr(self.recorder, table_name)
        self.insert(table_name, insert_query, func)

    def insert_imgt_inn_blastp(self) -> tuple:
        '''
        note: single one *.blastp5
        '''
        table_name = 'imgt_inn_blastp'
        insert_query = f"""INSERT INTO {table_name}
            (query_pro_id, inn_chain_id, ref_pro_id, pdb_chain_id,
                evalue, bit_score, identity)
            VALUES (%s, %s, %s, %s, %s, %s, %s)"""
        func = getattr(self.recorder, table_name)
        self.insert(table_name, insert_query, func, 2)

########################

    def insert_uniprot(self) -> tuple:
        table_name = 'uniprot'
        insert_query = f"""INSERT INTO {table_name}
            (uniprot_acc, name, organism, pro_len, pro_seq)
            VALUES (%s, %s, %s, %s, %s)"""
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

    def insert_uniprot_feature(self) -> tuple:
        table_name = 'uniprot_feature'
        insert_query = f"""INSERT INTO {table_name}
            (uniprot_acc, feature, start, end, seq, note, evidence)
            VALUES (%s, %s, %s, %s, %s, %s, %s)"""
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
    def insert_ab_aacdb(self):
        # AACDB
        table_name = 'ab_aacdb'
        insert_query = f"""INSERT INTO {table_name}
            (pdb_id, antibody_chains, antigen_chains,
                antibody, protein, drug, targets, method,
                organism, resolution, reference)
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
        """
        func = lambda rec: [rec,]
        self.insert(table_name, insert_query, func)

    def insert_meta_sabdab(self):
        # SAbDab
        table_name = 'meta_sabdab'
        insert_query = f"""INSERT INTO {table_name}
            (pdb_id, summary_file) VALUES (%s, %s)
        """
        func = lambda rec: [rec,]
        self.insert(table_name, insert_query, func)

    def insert_ab_sabdab(self):
        # SAbDab
        table_name = 'ab_sabdab'
        insert_query = f"""INSERT INTO {table_name}
            (pdb_id, model, heavy_chain, light_chain,
                antigen_chain, antigen_type, resolution)
            VALUES (%s, %s, %s, %s, %s, %s, %s)
        """
        func = lambda rec: [rec,]
        self.insert(table_name, insert_query, func)

    def insert_ab_abybank(self):
        # abYbank
        table_name = 'ab_abybank'
        insert_query = f"""INSERT INTO {table_name}
            (pdb_id, chain_type, chain_numbering)
            VALUES (%s, %s, %s)
        """
        func = lambda rec: [rec,]
        self.insert(table_name, insert_query, func)
    
    def insert_absd_source(self):
        table_name = 'absd_source'
        insert_query = f"""INSERT INTO {table_name}
            (seq_id, absd_id, header, dxref_id, vgene_family, source)
            VALUES (%s, %s, %s, %s, %s, %s)
        """
        func = getattr(self.recorder, table_name)
        self.insert(table_name, insert_query, func)

    def insert_absd_seq(self):
        table_name = 'absd_seq'
        insert_query = f"""INSERT INTO {table_name}
            (absd_id, specie, version, pro_len, pro_seq, record_id)
            VALUES (%s, %s, %s, %s, %s, %s)
        """
        func = lambda rec: [rec,]
        self.insert(table_name, insert_query, func)

    def insert_absd_faa(self):
        table_name = 'absd_faa'
        insert_query = f"""INSERT INTO {table_name}
            (seq_id, faa_file) VALUES (%s, %s)
        """
        func = lambda rec: [rec,]
        self.insert(table_name, insert_query, func)

    def insert_absd_igblastp(self):
        table_name = 'absd_igblastp'
        insert_query = f"""INSERT INTO {table_name}
            (seq_id, aln_file) VALUES (%s, %s)
        """
        func = lambda rec: [rec,]
        self.insert(table_name, insert_query, func)

    def insert_absd_align_imgt(self) -> tuple:
        table_name = 'absd_align_imgt'
        insert_query = f"""INSERT INTO {table_name}
            (seq_id, allele_name, gene_name, gene_family,
                chain_type, isotype, evalue, bit_score)
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s)
        """
        func = getattr(self.recorder, table_name)
        self.insert(table_name, insert_query, func)

    def insert_absd_align_vregion(self) -> tuple:
        table_name = 'absd_align_vregion'
        insert_query = f"""INSERT INTO {table_name}
            (seq_id, region_name, identity, start, end, seq)
            VALUES (%s, %s, %s, %s, %s, %s)
        """
        func = getattr(self.recorder, table_name)
        self.insert(table_name, insert_query, func)

#######################

    def insert_protein_seq(self):
        table_name = 'protein_seq'
        insert_query = f"""INSERT INTO {table_name}
            (pro_id, prefix, pro_len, pro_seq)
            VALUES (%s, %s, %s, %s)"""
        func = lambda rec: [rec,]
        self.insert(table_name, insert_query, func)

    def insert_protein_source(self):
        table_name = 'protein_source'
        insert_query = f"""INSERT INTO {table_name}
            (pro_id, source_name, source_id)
            VALUES (%s, %s, %s)"""
        func = lambda rec: [rec,]
        self.insert(table_name, insert_query, func)

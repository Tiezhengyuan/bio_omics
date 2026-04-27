'''
insert data into empty table
'''
import os
from .build_db import BuildDb
from .build_records import BuildRecords

# build SQL query string
TABLE_COLUMNS = {
    # meta
    'meta': ['pdb_group', 'host', 'ip_addr', 'script'],
    'meta_execute': ['table_name', 'pdb_group', 'host', 'ip_addr', 'script'],

    # 
    'pdb_protein': ['pdb_id', 'relative_fasta', 'pro_id', 'chain',
        'auth_chain', 'pro_len',  'pro_seq'],
    'meta_pdb': ['pdb_id', 'relative_path', 'pdb_group'],
    'meta_mmcif': ['pdb_id', 'relative_path', 'pdb_group'],
    'pdb_chain_type': ['pdb_id', 'chain_id', 'chain_type', 'total_res', 
        'non_water_res', 'aa_res', 'dna_res', 'rna_res', 'water_res',
        'ion_res', 'hetero_res', 'other_res'],
    'pdb_faa': ['pdb_id', 'model_no', 'relative_faa'],
    'meta_chain_faa': ['pdb_id', 'chain_id', 'relative_path'],
    'meta_chain_pdb': ['pdb_id', 'chain_id', 'chain_pdb'],
    'combo2': ['pdb_id', 'combo_id', 'chain_combo', 'chain_pdbs',
        'combo_pdb', 'combo_csv', 'combo_fasta'],
    'combo3': ['pdb_id', 'combo_id', 'chain_combo', 'chain_pdbs',
        'combo_pdb', 'combo_csv', 'combo_fasta'],
    'combo4': ['pdb_id', 'combo_id', 'chain_combo', 'chain_pdbs',
        'combo_pdb', 'combo_csv', 'combo_fasta'],
    'combo5': ['pdb_id', 'combo_id', 'chain_combo', 'chain_pdbs',
        'combo_pdb', 'combo_csv', 'combo_fasta'],

    'meta_combo2_pdb': ['pdb_id', 'model_no', 'combo_id',
        'chain_combo', 'bounded_pdb', 'chain_pdbs'],
    'meta_combo2_seq': ['pdb_id', 'model_no', 'combo_id',
        'chain_combo', 'relative_csv'],

    'label_chain': ['pdb_id', 'chain_id', 'first_chain_id',
        'chain_label', 'category'],
    'label_complex': ['pdb_id', 'model_no', 'combo_id', 'chain_combo', 
        'chain_label', 'combo', 'chain_ids', 'category'],

    #analyze
    'meta_bfactor': ['pdb_id', 'relative_pkl'],
    'combo2_distance': ['pdb_id', 'combo_id', 'relative_pkl',\
        'min_dist', 'min_res1', 'min_res2'],
    'pdb_distance': ['pdb_id', 'combo_id', 'relative_pkl'],
    'dist_amber14': ['pdb_id', 'combo_id', 'relative_pkl'],
    'combo2_sasa': ['pdb_id', 'combo_id', 'delta_total_sasa', 'relative_pkl'],
    'combo2_contacts': ['pdb_id', 'combo_id', 'binding_affinity',
        'dissociation_constant', 'relative_prodigy'],
    'pdb_sasa': ['pdb_id', 'combo_id', 'relative_pkl'],
    'equal_seq_rmsd': ['chain_id', 'equal_chain_id', 'rmsd',
        'chain_pdb', 'equal_chain_pdb',],

    # fix
    'fixed_pdb': ['pdb_id', 'relative_path', 'time0'],
    'fixed_chainid': ['pdb_id', 'model_no', 'fixed_model_no',
        'chain_no', 'fixed_chain_no'],
    'fixed_chain': ['pdb_id', 'chain_id', 'model_no',
        'chain_no', 'relative_pdb'],
    'fixed_combo2': ['pdb_id', 'model_no', 'combo_id',
        'chain_combo', 'relative_pdb'],
    'fixed_combo3': ['pdb_id', 'model_no', 'combo_id',
        'chain_combo', 'relative_pdb'],
    'pdb_chain_fixed': ['chain_id', 'relative_pdb', 'time0'],

    # simulate
    'relaxed_pdb_amber99sb': ['pdb_id', 'original_energy',
        'relaxed_energy', 'time0', 'relative_pdb'],
    'relaxed_pdb_amber14': ['pdb_id', 'original_energy',
        'relaxed_energy', 'time0', 'relative_pdb'],
    'relaxed_amber14_chain': ['pdb_id', 'chain_id', 'model_no',
        'chain_no', 'relative_pdb'],
    'relaxed_amber14_combo2': ['pdb_id', 'combo_id', 'model_no',
        'chain_combo', 'endpoint_prefix'],
    'relaxed_amber14_combo3': ['pdb_id', 'combo_id', 'model_no',
        'chain_combo', 'endpoint_prefix'],
    'chain_relaxed_amber14': ['chain_id', 'original_energy',
        'relaxed_energy', 'time0', 'endpoint'],

    'meta_chain_relaxed': ['pdb_id', 'chain_id', 'time0', 'fixed_pdb', 
        'relaxed_pdb', 'original_energy', 'relaxed_energy',],
    'combo2_relaxed': ['pdb_id', 'combo_id', 'time0', 'fixed_pdb', 
        'relaxed_pdb', 'original_energy', 'relaxed_energy',],
    'combo3_relaxed': ['pdb_id', 'combo_id', 'time0', 'fixed_pdb', 
        'relaxed_pdb', 'original_energy', 'relaxed_energy',],
    'combo4_relaxed': ['pdb_id', 'combo_id', 'time0', 'fixed_pdb', 
        'relaxed_pdb', 'original_energy', 'relaxed_energy',],
    'combo5_relaxed': ['pdb_id', 'combo_id', 'time0', 'fixed_pdb', 
        'relaxed_pdb', 'original_energy', 'relaxed_energy',],

    #align
    'meta_blastp5': ['pdb_id', 'model_no', 'relative_aln'],
    'meta_vdj': ['pdb_id', 'model_no', 'v_aln', 'd_aln', 'j_aln'],
    'align_dregion': ['pdb_id', 'chain_id', 'hit_id', 'organism',
        'evalue', 'bit_score', 'identities', 'query_len', 'align_len',
        'query_start', 'query_end', 'query_seq', 'match_seq',
        'sbjct_seq', 'sbjct_start', 'sbjct_end',],
    'align_jregion': ['pdb_id', 'chain_id', 'allele_name', 'specie',
        'gene_name', 'gene_family', 'chain_type', 'isotype',
        'evalue', 'bit_score', 'identities', 'query_len', 'align_len',
        'query_start', 'query_end', 'query_seq', 'match_seq',
        'sbjct_seq', 'sbjct_start', 'sbjct_end',],
    'ab_anarci': ['chain_id', 'domain_no', 'hmm_species',
        'chain_type', 'evalue', 'score', 'seqstart_index',
        'seqend_index', 'identity_species', 'v_gene', 'v_identity',
        'j_gene', 'j_identity', 'aligned_seq',],

    #predict
    'chain_alphafold': ['pdb_id', 'chain_id', 'endpoint', 'host', 'predicted_at',],
    'chain_plddt': ['chain_id', 'relative_pdb', 'ranking', 'predictor', 
        'avg_plddt', 'relative_pkl'],
    'chain_af_summary': ['chain_id', 'ranking', 'avg_plddt', 'avg_ptm',
        'max_pae', 'relative_json', 'relaxed_pdb', 'unrelaxed_pdb'],
    'chain_rmsd': ['pdb_id','chain_id', 'rmsd', 'chain_status', 
        'native_pdb', 'model_pdb'],
    'chain_tmalign': ['pdb_id','chain_id', 'tm1', 'tm2', 'rmsd',
        'chain_status', 'native_pdb', 'model_pdb'],
    'combo2_alphafold': ['pdb_id', 'combo_id', 'endpoint', 'host', 'predicted_at',],
    'combo3_alphafold': ['pdb_id', 'combo_id', 'endpoint', 'host', 'predicted_at',],


    # uniprot
    'uniprot_meta_alphafold': ['uniprot_acc', 'pdb_path'],
    'uniprot_alphafold_rmsd': ['chain_id', 'uniprot_acc', 'ca_rmsd', 
        'native_pdb', 'model_pdb'],
    'uniprot_rmsd': ['chain_id', 'native_pdb', 'uniprot_acc', 
        'model_pdb', 'ca_rmsd'],

    # imgt
    'imgt_gene': ['record_id', 'specie', 'isotype', 'chain_type',
        'gene_name', 'allele_name', 'region_name'],
    'imgt_geneseq': ['record_id', 'seq_type', 'seq'],
    'msa_vregion': ['chain_id', 'chain_types', 'query_fa',
        'aln_file', 'pickle_file'],

    # antibody
    'ab_aacdb': ['pdb_id', 'antibody_chains', 'antigen_chains',
        'antibody', 'protein', 'drug', 'targets', 'method',
        'organism', 'resolution', 'reference'],
    'meta_sabdab': ['pdb_id', 'summary_file'],
    'ab_sabdab': ['pdb_id', 'model', 'heavy_chain', 'light_chain',
        'antigen_chain', 'antigen_type', 'resolution'],
    'ab_abybank': ['pdb_id', 'chain_type', 'chain_numbering'],
    'absd_seq': ['absd_id', 'specie', 'version', 'pro_len', 
        'pro_seq', 'record_id'],
    'absd_faa': ['seq_id', 'faa_file'],
    'absd_igblastp': ['seq_id', 'aln_file'],

    #design
    'design_dimer': ['design_id', 'category', 'chain_labels', 'target_name',
        'target_id', 'binder_id', 'af_faa', 'target_pdb', 'binder_pdb'],
    'design_multimer': ['design_id', 'category', 'chain_labels', 'combo',
        'target_name', 'target_id', 'binder_id', 'af_faa', 'target_pdb',
        'binder_pdb'],
    'predict': ['predictor', 'design_id', 'host', 'status'],
}


class BuildComplex(BuildDb):

    def __init__(self, data:list=None, verbose:bool=False, chunk_size:int=None):
        '''
        args: data is iterable object of pdb_file, or record
        '''
        super().__init__(data, verbose, chunk_size)
        self.config['database'] = os.getenv('MYSQL_DB', 'COMPLEX2')
        self.recorder = BuildRecords(verbose)
    
    def insert_batch_records(self, table_name:str):
        cols = TABLE_COLUMNS[table_name]
        fmt = ['%s',] * len(cols)
        cols_str = ', '.join(cols)
        fmt_str= ', '.join(fmt)
        insert_query = f"""
            INSERT INTO {table_name} ({cols_str})
            VALUES ({fmt_str})
        """
        # load insertion
        func = lambda rec: [rec,]
        self.insert(table_name, insert_query, func)

    def insert_pdb(self):
        table_name = 'pdb'
        insert_query = f"""INSERT INTO {table_name}
            (pdb_id, pdb_url, release_date, name)
            VALUES (%s, %s, %s, %s)"""
        func = getattr(self.recorder, table_name)
        self.insert(table_name, insert_query, func)

    def insert_pdb_method(self):
        table_name = 'pdb_method'
        insert_query = f"""INSERT INTO {table_name}
            (pdb_id, structure_method, resolution)
            VALUES (%s, %s, %s)"""
        func = getattr(self.recorder, table_name)
        self.insert(table_name, insert_query, func)

    def insert_pdb_compound(self) -> tuple:
        table_name = 'pdb_compound'
        insert_query = f"""INSERT INTO {table_name} 
            (pdb_id, compound_no, chain_no, compound_text)
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
            (chain_id, pdb_id, model_no, chain_no, idx)
            VALUES (%s, %s, %s, %s, %s)"""
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
            (pdb_id, chain_id, pro_len, chain_seq)
            VALUES (%s, %s, %s, %s)"""
        func = getattr(self.recorder, table_name)
        return self.insert(table_name, insert_query, func)

    def insert_pdb_bfactor(self) -> tuple:
        table_name = 'pdb_bfactor'
        insert_query = f"""INSERT INTO {table_name} 
            (pdb_id, model_no, chain_no, avg_bfactor, bfactor)
            VALUES (%s, %s, %s, %s, %s)"""
        func = getattr(self.recorder, table_name)
        self.insert(table_name, insert_query, func)

    def insert_pdb_dssp(self) -> tuple:
        table_name = 'pdb_dssp'
        insert_query = f"""INSERT INTO {table_name} 
            (pdb_id, model_no, chain_no, dssp_seq)
            VALUES (%s, %s, %s, %s)"""
        func = getattr(self.recorder, table_name)
        self.insert(table_name, insert_query, func)

    def insert_align_uniprot(self) -> tuple:
        table_name = 'align_uniprot'
        insert_query = f"""INSERT INTO {table_name}
            (pdb_id, chain_id, uniprot_acc, evalue, bit_score, identity)
            VALUES (%s, %s, %s, %s, %s, %s)"""
        func = getattr(self.recorder, table_name)
        self.insert(table_name, insert_query, func, 2)

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
            (inn_number, chain_id, chain_description, 
                pro_len, chain_seq)
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
                gene_allele, seq, seq_cdr, seq_align,
                seq_align_cdr, seq_cdr_mask, seq_align_cdr_mask)
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
            (query_pro_id, inn_chain_id, ref_pro_id,
                pdb_chain_id, evalue, bit_score, identity)
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

    def insert_uniprot_alphafold_plddt(self) -> tuple:
        table_name = 'uniprot_alphafold_plddt'
        insert_query = f"""INSERT INTO {table_name}
            (uniprot_acc, model_no, chain_no, avg_plddt, pickle_file)
            VALUES (%s, %s, %s, %s, %s)"""
        func = getattr(self.recorder, table_name)
        self.insert(table_name, insert_query, func)

    def insert_align_vfrag(self) -> tuple:
        # FRs and CDRs in V-regions
        table_name = 'align_vfrag'
        insert_query = f"""INSERT INTO {table_name}
            (pdb_id, chain_id, identity, seq, fragment, seq_from, seq_to)
            VALUES (%s, %s, %s, %s, %s, %s, %s)
        """
        func = getattr(self.recorder, table_name)
        self.insert(table_name, insert_query, func)

    def insert_align_vregion(self) -> tuple:
        table_name = 'align_vregion'
        insert_query = f"""INSERT INTO {table_name}
            (pdb_id, chain_id, allele_name, specie, gene_name, gene_family,
                chain_type, isotype, evalue, bit_score, region_names)
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
        """
        func = getattr(self.recorder, table_name)
        self.insert(table_name, insert_query, func)

#####################
    
    def insert_absd_source(self):
        table_name = 'absd_source'
        insert_query = f"""INSERT INTO {table_name}
            (seq_id, absd_id, header, dxref_id, 
                vgene_family, source)
            VALUES (%s, %s, %s, %s, %s, %s)
        """
        func = getattr(self.recorder, table_name)
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

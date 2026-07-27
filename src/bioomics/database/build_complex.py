'''
insert data into empty table
'''
import os
from .insert_db import InsertDb
from .delete_db import DeleteDb
from .build_records import BuildRecords
'''
    
    #################
    #analyze
    'dist_amber14': ['pdb_id', 'combo_id', 'relative_pkl'],

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
'''

# build SQL query string
TABLE_COLUMNS = {
    # pdb
    'pdb': ['pdb_id', 'pdb_url', 'release_date', 'name',],
    'meta_pdb': ['pdb_id', 'relative_pdb', 'relative_cif', 'pdb_group'],
    'pdb_compound': ['pdb_id', 'compound_no', 'chain_no', 'compound_text',],
    'pdb_method': ['pdb_id', 'structure_method', 'resolution',],
    'pdb_keyword': ['pdb_id', 'keyword',],
    'pdb_protein': ['pdb_id', 'relative_fasta', 'pro_id', 'chain', 'auth_chain', 'pro_seq'],
    'meta_bfactor': ['pdb_id', 'relative_pkl'],
    'pdb_bfactor': ['pdb_id', 'model_no', 'chain_no', 'avg_bfactor', 'bfactor'],

    # pdb_chain
    'pdb_chainid': ['chain_id', 'pdb_id', 'model_no', 'chain_no', 'idx'],
    'pdb_chain_seq': ['pdb_id', 'chain_id', 'chain_seq'],
    'pdb_faa': ['pdb_id', 'model_no', 'relative_faa'],
    'chain_faa': ['pdb_id', 'chain_id', 'relative_path'],
    'pdb_pro_seq': ['pdb_id', 'chain_id', 'proseq_name', 'seq_id'],
    'pdb_chain_type': ['pdb_id', 'chain_id', 'chain_type', 'total_res', 
        'non_water_res', 'aa_res', 'dna_res', 'rna_res', 'water_res',
        'ion_res', 'hetero_res', 'other_res'],
    'chain_pdb': ['pdb_id', 'chain_id', 'chain_pdb'],
    'chain_dssp': ['pdb_id', 'chain_id', 'host', 'model_no',
        'chain_no', 'dssp_seq',],
   
    # uniprot
    'uniprot': ['uniprot_acc', 'pro_name', 'organism', 'pro_seq', 'creation_date', 'seq_date'],
    'uniprot_keyword': ['uniprot_acc', 'keyword'],
    'uniprot_dbxrefs': ['uniprot_acc', 'db_name', 'db_acc'],
    'uniprot_feature': ['uniprot_acc', 'feature', 'start', 'end', 
        'seq', 'note', 'evidence'],
    'uniprot_pro_seq': ['uniprot_acc', 'proseq_name', 'seq_id',],
    'uniprot_af': ['uniprot_acc', 'pdb_path'],
    # 'uniprot_alphafold_plddt': ['uniprot_acc', 'model_no', 'chain_no',
    #     'avg_plddt', 'pickle_file',],
    'uniprot_af_rmsd': ['chain_id', 'uniprot_acc', 'ca_rmsd', 
        'native_pdb', 'model_pdb'],
    # 'uniprot_rmsd': ['chain_id', 'native_pdb', 'uniprot_acc', 
    #     'model_pdb', 'ca_rmsd'],

    # label
    'meta_blastp5': ['pdb_id', 'model_no', 'relative_aln'],
    'align_uniprot': ['pdb_id', 'chain_id', 'uniprot_acc', 
        'evalue', 'bit_score', 'identity',],
    'meta_vdj': ['pdb_id', 'model_no', 'v_aln', 'd_aln', 'j_aln'],
    'align_vfrag': ['pdb_id', 'chain_id', 'identity', 'seq',
        'fragment', 'seq_from', 'seq_to'],
    'align_vregion':['pdb_id', 'chain_id', 'allele_name', 'specie',
        'gene_name', 'gene_family', 'chain_type', 'isotype',
        'evalue', 'bit_score', 'region_names'],
    'align_dregion': ['pdb_id', 'chain_id', 'hit_id', 'organism',
        'evalue', 'bit_score', 'identities', 'query_len', 'align_len',
        'query_start', 'query_end', 'query_seq', 'match_seq',
        'sbjct_seq', 'sbjct_start', 'sbjct_end',],
    'align_jregion': ['pdb_id', 'chain_id', 'allele_name', 'specie',
        'gene_name', 'gene_family', 'chain_type', 'isotype',
        'evalue', 'bit_score', 'identities', 'query_len', 'align_len',
        'query_start', 'query_end', 'query_seq', 'match_seq',
        'sbjct_seq', 'sbjct_start', 'sbjct_end',],
    'label_chain': ['pdb_id', 'chain_id', 'first_chain_id',
        'chain_label', 'category'],
    'label_complex': ['pdb_id', 'model_no', 'combo_id', 'chain_combo', 
        'chain_label', 'combo', 'chain_ids', 'category'],

    # dynamics.sh
    # combination
    'combo2': ['pdb_id', 'combo_id', 'chain_combo', 'chain_pdbs',
        'combo_pdb', 'combo_csv', 'combo_fasta'],
    'combo3': ['pdb_id', 'combo_id', 'chain_combo', 'chain_pdbs',
        'combo_pdb', 'combo_csv', 'combo_fasta'],
    'combo4': ['pdb_id', 'combo_id', 'chain_combo', 'chain_pdbs',
        'combo_pdb', 'combo_csv', 'combo_fasta'],
    'combo5': ['pdb_id', 'combo_id', 'chain_combo', 'chain_pdbs',
        'combo_pdb', 'combo_csv', 'combo_fasta'],
    # relax
    'chain_relaxed': ['pdb_id', 'chain_id', 'time0', 'fixed_pdb', 
        'relaxed_pdb', 'original_energy', 'relaxed_energy',],
    'combo2_relaxed': ['pdb_id', 'combo_id', 'time0', 'fixed_pdb', 
        'relaxed_pdb', 'original_energy', 'relaxed_energy',],
    'combo3_relaxed': ['pdb_id', 'combo_id', 'time0', 'fixed_pdb', 
        'relaxed_pdb', 'original_energy', 'relaxed_energy',],
    'combo4_relaxed': ['pdb_id', 'combo_id', 'time0', 'fixed_pdb', 
        'relaxed_pdb', 'original_energy', 'relaxed_energy',],
    'combo5_relaxed': ['pdb_id', 'combo_id', 'time0', 'fixed_pdb', 
        'relaxed_pdb', 'original_energy', 'relaxed_energy',],
    # contacts
    'equal_seq_rmsd': ['chain_id', 'equal_chain_id', 'rmsd',
        'chain_pdb', 'equal_chain_pdb',],
    'combo2_distance': ['pdb_id', 'combo_id', 'relative_pkl',\
        'min_dist', 'min_res1', 'min_res2'],
    'combo_sasa': ['pdb_id', 'combo_id', 'delta_total_sasa', 'relative_pkl'],
    'combo_contacts': ['pdb_id', 'combo_id', 'binding_affinity',
        'dissociation_constant', 'relative_prodigy'],

    #predict
    'chain_af': ['pdb_id', 'chain_id', 'endpoint', 'host', 'predicted_at',],
    'chain_af_summary': ['chain_id', 'ranking', 'avg_plddt', 'avg_ptm',
        'max_pae', 'relative_json', 'relaxed_pdb', 'unrelaxed_pdb'],
    'chain_af_plddt': ['chain_id', 'relative_pdb', 'ranking',
        'avg_plddt', 'relative_pkl'],
    'chain_rmsd': ['pdb_id','chain_id', 'rmsd', 'chain_status', 
        'native_pdb', 'model_pdb'],
    'chain_tmalign': ['pdb_id','chain_id', 'tm1', 'tm2', 'rmsd',
        'chain_status', 'native_pdb', 'model_pdb'],
    'combo2_af': ['pdb_id', 'combo_id', 'endpoint', 'host', 'predicted_at',],
    'combo3_af': ['pdb_id', 'combo_id', 'endpoint', 'host', 'predicted_at',],

    # meta
    'meta': ['pdb_group', 'host', 'ip_addr', 'script'],
    'meta_execute': ['table_name', 'pdb_group', 'host', 'ip_addr', 'script'],

    #design
    'design_chain': ['target_name', 'chain_id', 'parent_id', 'pdb_id',
        'chain_no', 'faa_file', 'pdb_file', 'binder_pdb',
        'chain_seq', 'epitope_pos'],
    'design_dimer': ['design_id', 'category', 'chain_labels', 'target_name',
        'target_id', 'target_chain_no', 'binder_id', 'binder_chain_no',
        'af_faa', 'target_pdb', 'binder_pdb', 'pdb_files'],
    'design_trimer': ['design_id', 'category', 'chain_labels', 'target_name',
        'target_id', 'target_chain_no', 'binder_id', 'binder_chain_no',
        'af_faa', 'target_pdb', 'binder_pdb', 'pdb_files'],
    
    'predict': ['target_name', 'predictor', 'design_id', 'host', 'status'],
    'predict_summary': ['summary_id', 'predict_id', 'ranking', 'avg_plddt',
        'avg_ptm', 'max_pae', 'relative_json', 'relaxed_pdb', 'unrelaxed_pdb'],
    'unrelaxed_plddt': ['summary_id', 'relative_pdb', 'relative_pkl', 
        'chain_no', 'avg_plddt',],
    'relaxed_plddt': ['summary_id', 'relative_pdb', 'relative_pkl', 
        'chain_no', 'avg_plddt',],
    'unrelaxed_contacts': ['summary_id', 'relative_pdb', 'binding_affinity',
        'dissociation_constant', 'relative_prodigy'],
    'relaxed_contacts': ['summary_id', 'relative_pdb', 'binding_affinity',
        'dissociation_constant', 'relative_prodigy'],
    'relaxed_dist': ['summary_id', 'relative_pdb', 'chain_combo',
        'min_dist', 'min_res1', 'min_res2', 'relative_pkl'],
    'unrelaxed_dist': ['summary_id', 'relative_pdb', 'chain_combo',
        'min_dist', 'min_res1', 'min_res2', 'relative_pkl'],
    'unrelaxed_sasa': ['summary_id', 'delta_total_sasa', 'relative_pkl', 'relative_pdb',],
    'relaxed_sasa': ['summary_id', 'delta_total_sasa', 'relative_pkl', 'relative_pdb',],
    'unrelaxed_dist_cdr': ['summary_id', 'binder_id', 'fragment', 
        'res_no', 'res', 'min_dist',],
    'relaxed_dist_cdr': ['summary_id', 'binder_id', 'fragment', 
        'res_no', 'res', 'min_dist',],


    'docking': ['target_name', 'design_id', 'host', 'status'],
    


    'ab_anarci': ['chain_id', 'domain_no', 'hmm_species',
        'chain_type', 'evalue', 'score', 'seqstart_index',
        'seqend_index', 'identity_species', 'v_gene', 'v_identity',
        'j_gene', 'j_identity', 'aligned_seq',],
}


class BuildComplex(InsertDb, DeleteDb):

    def __init__(self, verbose:bool=False, chunk_size:int=None):
        '''
        args: data is iterable object of pdb_file, or record
        '''
        InsertDb.__init__(self, verbose, chunk_size)
        DeleteDb.__init__(self, verbose)
        self.config['database'] = os.getenv('MYSQL_DB', 'COMPLEX2')
        self.recorder = BuildRecords(verbose)
    
    def insert_batch_records(self, data, table_name:str, cols=None):
        if cols is None:
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
        if self.verbose:
            print(f"\nBegin to insert records into {table_name}...", end='')
        return self.func_insert(data, insert_query, func)

   
    def insert_pdb_chain(self, data):
        table_name = 'pdb_chain'
        insert_query = f"""INSERT INTO {table_name} 
            (chain_id, pdb_id, model_no, chain_no, chain_seq)
            VALUES (%s, %s, %s, %s, %s)"""
        func = getattr(self.recorder, table_name)
        return self.func_insert(data, insert_query, func)

    def insert_proseq(self, rows:list, source=None) -> tuple:
        # separate sequences by length
        seq400, seq1000, seqlong = [], [], []
        for row in rows:
            seq = row['seq']
            rec = (source, seq)
            # duplicates should be excluded
            if len(seq) <= 400:
                if rec not in seq400:
                    seq400.append(rec)
            elif 400 < len(seq) <= 1000:
                if rec not in seq1000:
                    seq1000.append(rec)
            else:
                if rec not in seqlong:
                    seqlong.append(rec)
        #
        cols = ['source', 'seq']
        s1, f1 = self.insert_records(seq400, 'proseq_400', cols)
        s2, f2 = self.insert_records(seq1000, 'proseq_1000', cols)
        s3, f3 = self.insert_records(seqlong, 'proseq_long', cols)
        return s1 + s2 + s3, f1 + f2 + f3
    
    def put_design_chain(self, data):
        table_name = 'design_chain'

        # delete old rows
        target_name = data['target_name']
        chain_id = data['chain_id']
        query = f"""
            DELETE FROM {table_name}
            WHERE target_name = '{target_name}'
                AND chain_id = '{chain_id}'
        """
        self.execute(query)

        # insert new record
        cols = TABLE_COLUMNS[table_name]
        rec = [data.get(i) for i in cols]
        if self.verbose:
            print(rec)
        self.insert_batch_records([rec,], table_name)
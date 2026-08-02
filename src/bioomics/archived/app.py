#! /usr/bin/python
'''
example:
    python src/collect/distance.py
'''
import os
import sys
import subprocess

pool = {
    #pdb
    'download_pdb': 'pdb.download_pdb',
    'pdb': 'pdb.pdb',
    'meta_pdb': 'pdb.meta_pdb',
    'pdb_compound': 'pdb.pdb_compound',
    'pdb_method': 'pdb.pdb_method',
    'pdb_keyword': 'pdb.pdb_keyword',
    'pdb_chainid': 'pdb.pdb_chainid',
    'pdb_protein': 'pdb.pdb_protein',
    # pdb-analysis
    'meta_bfactor': 'pdb.meta_bfactor',
    'pdb_bfactor': 'pdb.pdb_bfactor',
    'chain_dssp': 'pdb.chain_dssp',

    # pdb-chain
    'chain_seq': 'pdb.chain_seq',
    'unique_chain_seq': 'pdb.unique_chain_seq',
    'pdb_faa': 'pdb.pdb_faa',
    'chain_faa': 'pdb.chain_faa',
    'pdb_chain_type': 'pdb.pdb_chain_type',
    'pdb_pro_seq': 'pdb.pdb_pro_seq',
    'chain_pdb': 'pdb.chain_pdb',


    # uniprot
    'uniprot': 'uniprot.uniprot',
    'uniprot_keyword': 'uniprot.uniprot_keyword',
    'uniprot_dbxrefs': 'uniprot.uniprot_dbxrefs',
    'uniprot_feature': 'uniprot.uniprot_feature',
    'uniprot_pro_seq': 'uniprot.uniprot_pro_seq',
    'uniprot_af': 'uniprot.uniprot_af',
    'uniprot_af_rmsd': 'uniprot.uniprot_af_rmsd',

    #dynamics
    'chain_relaxed': 'dynamics.chain_relaxed',
    'equal_seq_rmsd': 'dynamics.equal_seq_rmsd',
    'combo2': 'dynamics.combo2',
    'combo2_relaxed': 'dynamics.combo2_relaxed',
    'combo2_distance': 'dynamics.combo2_distance',
    'combo3': 'dynamics.combo3',
    'combo3_relaxed': 'dynamics.combo3_relaxed',
    'combo4': 'dynamics.combo4',
    'combo4_relaxed': 'dynamics.combo4_relaxed',
    'combo5': 'dynamics.combo5',
    'combo5_relaxed': 'dynamics.combo5_relaxed',

    'combo_sasa': 'dynamics.combo_sasa',
    'combo_contacts': 'dynamics.combo_contacts',

    # align
    'meta_blastp5': 'align.meta_blastp5',
    'align_uniprot': 'align.align_uniprot',
    'meta_vdj': 'align.meta_vdj',
    'align_vfrag': 'align.align_vfrag',
    'align_vregion': 'align.align_vregion',
    'align_dregion': 'align.align_dregion',
    'align_jregion': 'align.align_jregion',
    'label_antibody': 'align.label_antibody',
    # 'anarci': 'immune.anarci',

    # predict
    'chain_af': 'predict.chain_af',
    'chain_af_plddt': 'predict.chain_af_plddt',
    'chain_af_summary': 'predict.chain_af_summary',
    'chain_rmsd': 'predict.chain_rmsd',
    'chain_tmalign': 'predict.chain_tmalign',
    'combo2_af': 'predict.combo2_af',
    'combo3_af': 'predict.combo3_af',
    
    # design
    'predict_summary': 'design.predict_summary',
    'predict_plddt': 'design.predict_plddt',
    'predict_contacts': 'design.predict_contacts',
    'predict_dist': 'design.predict_dist',
    'predict_sasa': 'design.predict_sasa',
    'dist_cdr': 'design.dist_cdr',

    # collect data
    'collect_distance': 'collect.collect_distance',
    'collect_ab_distance': 'collect.collect_ab_distance',

    #other
    'clean_predict': 'clean.clean_predict',
    'move': 'clean.move',
    'test': 'test',
}

# def main():
#     # build cmd
#     project_dir = os.path.dirname(__file__)
#     name = sys.argv[1] if len(sys.argv) > 1 else None


#     cmd = ['python', '-m',]
#     path_prefix = 'src.pipelines.'
#     if name in pool:
#         cmd.append(path_prefix + pool[name])
#     else:
#         print("\nError: The argument <module_name> is not defined.\n")
#         print("Usage: python app.py <module_name> [pdb_group]\n")
#         sys.exit(1)

#     # pass arguments
#     args = sys.argv[2:] if len(sys.argv) > 2 else []
#     cmd += args
#     print(f"Run python script: {cmd}")

#     # run
#     result = subprocess.run(cmd)


def main():
    print('ccc')
if __name__ == "__main__":
    main()
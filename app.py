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
    'pdb_dssp': 'pdb.pdb_dssp',

    # pdb-seq
    'chain_seq': 'pdb.pdb_chain_seq',
    'unique_chain_seq': 'pdb.unique_chain_seq',
    'pdb_faa': 'pdb.pdb_faa',
    'meta_chain_faa': 'pdb.meta_chain_faa',
    'pdb_chain_type': 'pdb.pdb_chain_type',
    'protein_pdb': 'pdb.protein_pdb',

    #uniprot
    'blastp_uniprot': 'protein.blastp_uniprot',
    'align_uniprot': 'protein.align_uniprot',

    #IMGT
    'meta_vdj': 'immune.meta_vdj',
    'align_vfrag': 'immune.align_vfrag',
    'align_vregion': 'immune.align_vregion',
    'align_dregion': 'immune.align_dregion',
    'align_jregion': 'immune.align_jregion',
    'pdb_complex_antibody': 'immune.pdb_complex_antibody',
    'anarci': 'immune.anarci',

    #structures
    'meta_chain_pdb': 'pdb.meta_chain_pdb',
    'combo2': 'pdb.combo2',
    'combo3': 'pdb.combo3',
    'combo4': 'pdb.combo4',
    'combo5': 'pdb.combo5',

    #dynamics
    'meta_chain_relaxed': 'dynamics.meta_chain_relaxed',
    'combo2_relaxed': 'dynamics.combo2_relaxed',
    'combo2_distance': 'dynamics.combo2_distance',
    'combo2_sasa': 'dynamics.combo2_sasa',
    'equal_seq_rmsd': 'dynamics.equal_seq_rmsd',
    'combo2_contacts': 'dynamics.combo2_contacts',
    'combo3_relaxed': 'dynamics.combo3_relaxed',
    'combo4_relaxed': 'dynamics.combo4_relaxed',
    'combo5_relaxed': 'dynamics.combo5_relaxed',

    # collect data
    'collect_distance': 'collect.collect_distance',
    'collect_ab_distance': 'collect.collect_ab_distance',

    # predict
    'chain_af': 'predict.chain_af',
    'chain_af_plddt': 'predict.chain_af_plddt',
    'chain_af_summary': 'predict.chain_af_summary',
    'chain_rmsd': 'predict.chain_rmsd',
    'chain_tmalign': 'predict.chain_tmalign',
    'combo2_af': 'predict.combo2_af',
    'combo3_af': 'predict.combo3_af',

    #other
    'clean_predict': 'clean.clean_predict',
    'move': 'clean.move',
    'test': 'test',
}

# build cmd
project_dir = os.path.dirname(__file__)
name = sys.argv[1] if len(sys.argv) > 1 else None


cmd = ['python', '-m',]
path_prefix = 'src.pipelines.'
if name in pool:
    cmd.append(path_prefix + pool[name])
else:
    print("\nError: The argument <module_name> is not defined.\n")
    print("Usage: python app.py <module_name> [pdb_group]\n")
    sys.exit(1)

# pass arguments
args = sys.argv[2:] if len(sys.argv) > 2 else []
cmd += args
print(f"Run python script: {cmd}")

# run
result = subprocess.run(cmd)

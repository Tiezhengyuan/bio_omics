#! /usr/bin/python
'''
example:
    python src/collect/distance.py
'''
import sys
import subprocess

# build cmd
name = sys.argv[1] if len(sys.argv) > 1 else None
pdb_group = sys.argv[2] if len(sys.argv) > 2 else None
pool = {
    #pdb
    'download_pdb': 'pdb.download_pdb',
    'pdb': 'pdb.pdb',
    'meta_pdb': 'pdb.meta_pdb',
    'pdb_compound': 'pdb.pdb_compound',
    'pdb_method': 'pdb.pdb_method',
    'pdb_keyword': 'pdb.pdb_keyword',
    'pdb_chainid': 'pdb.pdb_chainid',
    'chain_seq': 'pdb.pdb_chain_seq',
    'unique_chain_seq': 'protein.unique_chain_seq',
    'pdb_faa': 'pdb.pdb_faa',
    'meta_chain_faa': 'pdb.meta_chain_faa',
    'meta_chain_pdb': 'pdb.meta_chain_pdb',

    #immune
    'meta_vdj': 'immune.meta_vdj',
    'align_vfrag': 'immune.align_vfrag',
    'align_dregion': 'immune.align_dregion',
    'align_jregion': 'immune.align_jregion',
    'anarci': 'immune.anarci',
    'pdb_complex': 'pdb.pdb_complex',

    #dynamics
    'meta_combo2': 'pdb.meta_combo2',
    'meta_combo3': 'pdb.meta_combo3',
    'combo2_distance': 'dynamics.combo2_distance', 
    'equal_seq_rmsd': 'dynamics.equal_seq_rmsd',
    'combo2_relaxed': 'dynamics.combo2_relaxed',
    'combo2_contacts': 'dynamics.combo2_contacts',
    'combo3_relaxed': 'dynamics.combo3_relaxed',

    # collect data
    'collect_distance': 'collect.collect_distance',
    'collect_ab_distance': 'collect.collect_ab_distance',

    # predict
    'chain_alpha': 'predict.chain_alphafold',
    'chain_alpha_plddt': 'predict.chain_alpha_plddt',
    'chain_afsum': 'predict.chain_afsum',
    'chain_rmsd': 'predict.chain_rmsd',
    'chain_tmalign': 'predict.chain_tmalign',
    'combo2_alpha': 'predict.combo2_alpha',
    'combo3_alpha': 'predict.combo3_alpha',

    #other
    'clean_predict': 'clean.clean_predict',
    'test': 'test',
}

cmd = ['python', '-m',]
path_prefix = 'src.pipelines.'
if name in pool:
    cmd.append(path_prefix + pool[name])
    if pdb_group is not None:
        cmd.append(pdb_group)
else:
    print("\nError: The argument <module_name> is not defined.\n")
    print("Usage: python app.py <module_name> [pdb_group]\n")
    sys.exit(1)

print(f"Run python script: {cmd}")
result = subprocess.run(cmd)

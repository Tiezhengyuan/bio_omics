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
    #immune
    'anarci': 'immune.anarci',

    #dynamics
    'combo2_distance': 'dynamics.combo2_distance', 
    'equal_seq_rmsd': 'dynamics.equal_seq_rmsd',
    'combo2_relaxed': 'dynamics.combo2_relaxed',
    'combo2_contacts': 'dynamics.combo2_contacts',

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

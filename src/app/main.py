#! /usr/bin/python
import os
import sys
import subprocess

pool = {
    #pdb
    'download_pdb': 'pdb.download_pdb',
    'download_cif': 'pdb.download_cif',
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

    # pdb-chain
    'pdb_chain_seq': 'pdb.pdb_chain_seq',
    'unique_chain_seq': 'pdb.unique_chain_seq',
    'pdb_faa': 'pdb.pdb_faa',
    'pdb_chain_type': 'pdb.pdb_chain_type',
    'pdb_proseq': 'pdb.pdb_proseq',
    'chain_pdb': 'pdb.chain_pdb',


    # uniprot
    'download_uniprot': 'uniprot.download_uniprot',
    'uniprot': 'uniprot.uniprot',
    'uniprot_keyword': 'uniprot.uniprot_keyword',
    'uniprot_dbxrefs': 'uniprot.uniprot_dbxrefs',
    'uniprot_feature': 'uniprot.uniprot_feature',
    'uniprot_proseq': 'uniprot.uniprot_proseq',
    'uniprot_align_db': 'uniprot.uniprot_align_db',
    'uniprot_af': 'uniprot.uniprot_af',
    'uniprot_af_rmsd': 'uniprot.uniprot_af_rmsd',

    #dynamics
    'equal_seq_rmsd': 'dynamics.equal_seq_rmsd',
    'combo2': 'dynamics.combo2',
    'combo3': 'dynamics.combo3',
    'combo4': 'dynamics.combo4',
    'combo5': 'dynamics.combo5',

    # align
    'meta_blastp5': 'align.meta_blastp5',
    'align_uniprot': 'align.align_uniprot',
    # align_vdj/vfrag/vregion/dregion/jregion, label_antibody moved to antibody_omics
    # 'anarci': 'immune.anarci',

    # predict
    'af_params': 'predict.af_params',
    'af_msa': 'predict.af_msa',
    'af_pred': 'predict.af_pred',
    'af_summary': 'predict.af_summary',
    'af_plddt': 'predict.af_plddt',
    'af_rmsd': 'predict.af_rmsd',
    'af_tmalign': 'predict.af_tmalign',

    
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
    'clean_empty': 'clean.clean_empty',
    'clean_predict': 'clean.clean_predict',
    'move': 'clean.move',
    'test': 'test',
}

def main():
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

if __name__ == "__main__":
    main()
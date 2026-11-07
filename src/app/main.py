#! /usr/bin/python
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

    # dynamics
    'equal_seq_rmsd': 'pdb.equal_seq_rmsd',

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
    'meta_blastp5': 'uniprot.meta_blastp5',
    'align_uniprot': 'uniprot.align_uniprot',

    # dynamics: combo2..5 moved to structure_profiling

    # align
    # meta_blastp5, align_uniprot moved to uniprot
    # align_vdj/vfrag/vregion/dregion/jregion, label_antibody moved to antibody_omics
    # 'anarci': 'immune.anarci',

    # predict: af_params, af_msa, af_pred, af_summary, af_plddt, af_rmsd,
    # af_tmalign moved to structure_profiling:
    #   see ~/bio/structure_profiling/bin/6_predict.sh

    
    # design: predict_summary, predict_plddt, predict_contacts, dist_cdr
    # predict_dist / predict_sasa moved to denovo_design:
    #   see ~/bio/denovo_design/pipelines/
    # collect_distance moved to antibody_omics:
    #   see ~/bio/antibody_omics/src/pipelines/collect/

    # collect_ab_distance moved to complex_contacts:
    #   see ~/bio/complex_contacts/src/pipelines/collect_ab_distance.py

    #other
    'clean_empty': 'clean.clean_empty',
    'clean_predict': 'clean.clean_predict',
    'move': 'clean.move',
    'test': 'test',
}

def main():
    # build cmd
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
    subprocess.run(cmd)

if __name__ == "__main__":
    main()
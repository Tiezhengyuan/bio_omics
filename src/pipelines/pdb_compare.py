'''
example: 
    - parallel -j4 python /home/yuan/bio/bio_omics/src/pipelines/pdb_compare.py \
        ::: /home/yuan/data/pdb/_list/pdb_*.json
    - python /home/yuan/bio/bio_omics/src/pipelines/pdb_compare.py \
        /home/yuan/data/pdb/_list/pdb_46.json
functions:
    - similarity between model and sample structure
'''
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import ProcessJson, QueryComplex, AnalyzeModel, BuildComplex


def parse_args(args_pool):
    rows = QueryComplex().parse_modelling()
    res = []
    for args in args_pool:
        for row in rows:
            if args['pdb_id'] == row['pdb_id']:
                args.update(row)
                res.append(args)
                break
    print(len(res))
    return res

def uniprot_alphafold(args):
    # retrieve chain pdb in 
    # retrieve uniprot pdb predicted by alphafold
    pdb_id = args['pdb_id']
    chain_id = args['chain_id']
    chain_no = chain_id.split('_')[-1]
    native_pdb = os.path.join(args['outdir'], \
        pdb_id[:2], pdb_id, 'chains', f'{pdb_id}_{chain_no}.pdb')
    if not os.path.isfile(native_pdb):
        print(args, native_pdb)
        return 'null'
    model_pdb = args['pdb_path']
    if not os.path.isfile(model_pdb):
        print(args)
        return 'null'

    # calculate RMSD
    try:
        rmsd = AnalyzeModel.rmsd(native_pdb, model_pdb)
        if rmsd:
            rec = (chain_id, native_pdb, args['uniprot_acc'], model_pdb, rmsd)
            BuildComplex(verbose=args['verbose']).insert_uniprot_rmsd(rec)
            return 'succeed'
    except Exception as e:
        pass
    return 'fail'


if __name__ == "__main__":
    status = ['succeed', 'skip', 'null', 'na', 'fail']
    meta = {
        'uniprot_alphafold': {i:0 for i in status},
    }
    json_file = sys.argv[1]
    args_pool = ProcessJson.load_json(json_file)
    args_pool = parse_args(args_pool)

    for args in args_pool:
        args['verbose'] = True
        args['overwrite'] = False
        status = uniprot_alphafold(args)
        meta['uniprot_alphafold'][status] += 1
    print(meta)




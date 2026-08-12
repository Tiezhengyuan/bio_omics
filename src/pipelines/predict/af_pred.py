'''
requirements:
    - tables af_msa
example: 
    - parallel -j8 bomics af_pred ::: {0..241}
    - bomics af_pred 0
functions:
    - predict structure
    - retrieve summary of alphafold prediction to table af_summary
'''

from src.helper import *
from bioomics import ProcessJson, Predictor
from .func_predict import *
from .func_summary import *



            
def parse_data(row):
    res = {}
    for path in row['outdir'].iterdir():
        if 'rank_' in path.name:
            ranking = int(re.findall(r'rank_(\d+)_', path.name)[0])
            if ranking not in res:
                res[ranking] = {'ranking': ranking,}
            if '_scores_' in path.name and path.suffix == '.json':
                data = ProcessJson.load_json(path)
                res[ranking]['avg_plddt'] = np.mean(data['plddt'])
                res[ranking]['avg_ptm'] = np.mean(data['ptm'])
                res[ranking]['max_pae'] = data['max_pae']
                res[ranking]['json_file_name'] = path.name
            elif '_relaxed_' in path.name and path.suffix == '.pdb':
                res[ranking]['relaxed_pdb_name'] = path.name
            elif '_unrelaxed_' in path.name and path.suffix == '.pdb':
                res[ranking]['unrelaxed_pdb_name'] = path.name
    return res

def build(row, pdb_names, meta):
    records = []
    ranked = parse_data(row)
    for ranking, log in ranked.items():
        relaxed_pdb_name = log.get('relaxed_pdb_name')
        if relaxed_pdb_name not in pdb_names:
            rec = (
                row['input_id'],
                ranking,
                log.get('avg_plddt'),
                log.get('avg_ptm'),
                log.get('max_pae'),
                log.get('json_file_name'),
                relaxed_pdb_name, 
                log.get('unrelaxed_pdb_name'),
            )
            records.append(rec)
    meta['records'] += len(ranked)
    return records



if __name__ == "__main__":
    params.update({
        'tmp_dir': os.getenv('tmp_dir'),
        'predictor': 'AF2',
        'num_recycle': '3',
        'relax': True,
        'params_table_name': 'pdb_pred_params',
        'msa_table_name': 'af_msa',
        'table_name': 'af_pred',
    })
    
    # a3m files
    rows = retrieve_msa(params, meta)
    
    # detect files in af_pred table
    relaxed_pdb_names = pred_summary(params)

    bc = BuildComplex(params['verbose'])
    for row in rows:
        run_colabfold_batch(row, params, meta)
        
        # parse results
        records = build(row, relaxed_pdb_names, meta)
        if records:
            print(f"Try to insert {row['input_id']} into {params['table_name']}")
            bc.insert_batch_records(records, params['table_name'])
    footer(meta)
    
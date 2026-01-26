'''
requirements:
    - alphafold prediciton model in PDB
    - table uniprot_meta_alphafold
example: 
    - python src/pipelines/uniprot_alphafold_plddt.py
functions:
    - calculate pLDDT by uniprot protein
    - put pLDDT to table uniprot_alphafold_plddt
'''

from helper import *
from bioomics import ProcessPickle, CalculatePlddt

def retrieve_rows(meta):
    query = f"""
        SELECT uniprot_acc, pdb_path
        FROM uniprot_meta_alphafold
        WHERE uniprot_acc NOT IN (
                SELECT uniprot_acc
                FROM uniprot_alphafold_plddt
            );
    """
    # each record must have pdb_id and pdb_file
    rows = QueryComplex().list_data(query)
    meta['rows'] = len(rows)
    print(f"{len(rows)} PDB files ready...")
    return rows

def cal_plddt(rows, args, meta):
    outdir = os.path.join(args['data_dir'], 'uniprot_alphafold_plddt')
    Dir(outdir).init_dir()

    for row in rows:
        data = None
        uniprot_acc = row['uniprot_acc']
        pfile = os.path.join(outdir, f'{uniprot_acc}.pkl')
        if args['overwrite'] or (not os.path.isfile(pfile)):
            pdb_file = row['pdb_path']
            # load structure before calculation
            cp = CalculatePlddt(pdb_file, outdir)
            cp.load_structure()
            if args['verbose']:
                print(f"Try to calclate p given {pdb_file}...")
            data = cp.complex_ca_plddt()
            if data.get('pLDDT'):
                if args['verbose']:
                    print(f"Save pLDDT results to {pfile}")
                status = ProcessPickle(pfile).save_pickle(data)
                meta[status] += 1
            else:
                meta['null'] += 1
        else:
            meta['skip'] += 1
        # prepare data for insertion
        if data and data.get('pLDDT'):
            data['pfile'] = pfile
            data['uniprot_acc'] = uniprot_acc
            yield data

if __name__ == "__main__":
    args = {
        'verbose': True,
        'overwrite': False,
        'data_dir': os.getenv('data_dir'),
    }
    # retrieve pdb_ids
    rows = retrieve_rows(meta)

    # calculate pLDDT
    plddt_data_iter = cal_plddt(rows, args, meta)

    print("Try to insert pLDDT data into table")
    # insert data
    bc = BuildComplex(plddt_data_iter, args['verbose'])
    bc.insert_uniprot_alphafold_plddt()

    print(meta)

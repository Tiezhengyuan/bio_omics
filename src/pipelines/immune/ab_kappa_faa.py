'''
example:
    - python src/pipelines/ab_kappa_faa.py
functions:
    - collect antibody chain sequences and export to collect/
    - do multiple alignment
'''

from helper import *
from bioomics import AlignSeq, ProcessFasta

def pull_data(args, meta):
    query = '''
        SELECT chain_id AS id, chain_seq AS seq
        FROM view_antibody
        WHERE chain_type = 'K'
    '''
    df = QueryComplex(args['verbose']).list_data(query, True)
    meta['num_records'] = len(df)
    return df

if __name__ == "__main__":
    args = {
        'verbose': True,
        'data_dir': os.getenv('data_dir'),
        'muscle_bin': os.getenv("muscle_bin"),
    }
    args['collect_dir'] = os.path.join(args['data_dir'], 'collect')
    Dir(args['collect_dir']).init_dir()

    # chain seq
    df = pull_data(args, meta)

    # prepare fasta
    name = 'antibody_kappa'
    fa_file = os.path.join(args['collect_dir'], f'{name}.faa')
    ProcessFasta(fa_file).unique_sequences(df)

    # do MSA
    aln_file = os.path.join(args['collect_dir'], f'{name}.muscle')
    AlignSeq(args['muscle_bin'], args['verbose']).muscle(fa_file, aln_file, True)

    print(meta)


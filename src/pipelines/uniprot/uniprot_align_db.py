'''
requirements:
    - uniprot_sprot.dat is downloaded
example:
    - bomics uniprot_align_db
function:
    - download data from ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/
'''

from src.helper import *
from bioomics import ProcessBlast
    
def blastp_db_uniprot_sprot(params):
    params['blastp_db_uniprot_sprot'] = params['rawdata_uniprot_dir'] / 'db'
    params['blastp_db_uniprot_sprot'].mkdir(parents=True, exist_ok=True)
    pb = ProcessBlast(params['blastp_db_uniprot_sprot'], params['verbose'])
    
    fa_file = params['rawdata_uniprot_dir'] / "uniprot_sprot.fasta.gz"
    db_path = pb.build_blastp_db(fa_file)
    print('Blastp database:', db_path)
   
if __name__ == "__main__":
    params.update({
        'overwrite': False,
    })
    
    blastp_db_uniprot_sprot(params)    

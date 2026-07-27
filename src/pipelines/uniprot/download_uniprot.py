'''
requirements:
    - uniprot_sprot.dat is downloaded
example:
    - bomics download_uniprot
function:
    - download data from ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/
'''

from src.helper import *
from bioomics import DownloadUniprot
       
   
if __name__ == "__main__":
    params.update({
        'overwrite': False,
        'run_gunzip': False,
    })

    du =DownloadUniprot(params['rawdata_dir'], params['overwrite'], params['run_gunzip'])
    
    print(f"Try to download UniprotSprot...", end='')
    res = du.uniprot_sprot()
    print(f"Done! local paths= {res}")
    
    # print(f"Try to download UniprotTrembl...", end='')
    # res = du.uniprot_trembl()
    # print(f"Done! local paths= {res}")
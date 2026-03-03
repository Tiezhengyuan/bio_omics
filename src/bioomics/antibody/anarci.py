'''
ANARCI -i ~/bio/ANARCI/Example_scripts_and_sequences/12e8.fasta \
    --outfile ~/bio/ANARCI/Example_scripts_and_sequences/12e8 --csv

'''
import os
import pandas as pd
import subprocess

class Anarci:

    def __init__(self, outdir:str=None, verbose:bool=None):
        self.outdir = outdir
        self.verbose = verbose if verbose else False
    
    def run_anarci(self, query_fa:str, file_name:str=None):
        cmd = ['ANARCI', '-i', query_fa,]
        if os.path.isdir(self.outdir):
            if not file_name:
                file_name = 'result'
            outfile = os.path.join(self.outdir, file_name)
            cmd += ['--outfile', outfile, '--csv',]
        try:
            cmd = ' '.join(cmd)
            print(cmd)
            result = subprocess.run(cmd, text=True, check=True,
                capture_output=True, shell=True)
            if self.verbose:
                print("Stdout:", result.stdout)
                print("Stderr:", result.stderr)
                print("Exit code:", result.returncode)
            return 'succeed'
        except Exception as e:
            if self.verbose:
                print(query_fa, e)
        return 'fail'

    def parse_numbering(self, csv_files:list):
        '''
        cols index 0:12
            'Id', 'domain_no', 'hmm_species', 'chain_type', 'e-value',
            'score', 'seqstart_index', 'seqend_index', 'identity_species',
            'v_gene', 'v_identity', 'j_gene', 'j_identity',
        '''
        res = []
        for infile in csv_files:
            if os.path.isfile(infile) and infile.endswith('.csv'):
                df = pd.read_csv(infile)
                for _, row in df.iterrows():
                    rec = row[:13].to_dict()
                    rec['aligned_seq'] = row[13:].str.cat(sep='')
                    res.append(rec)
        return res

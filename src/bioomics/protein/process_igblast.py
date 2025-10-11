'''
install blast: sudo apt install ncbi-blast+
'''
import os
import sys
import subprocess
from biosequtils import Dir

class ProcessIgblast:

    def __init__(self, outdir:str=None, verbose:bool=False):
        self.outdir = outdir
        Dir(outdir).init_dir()
        self.verbose = verbose

    def build_prot_db(self, fa_file:str, bin_dir:str=None):
        outdir = self.outdir if os.path.isdir(self.outdir) else os.path.dirname(fa_file)
        prefix = os.path.splitext(os.path.basename(fa_file))[0]
        outprefix = os.path.join(outdir, prefix)
        
        if bin_dir is None:
            bin_dir = '/home/yuan/data/ncbi-igblast-1.22.0/bin'
        if bin_dir not in sys.path:
            sys.path.append(bin_dir)
        exe = os.path.join(bin_dir, 'makeblastdb') if bin_dir else 'makeblastdb'

        cmd = [ exe, '-parse_seqids', '-dbtype', 'prot', \
            '-in', fa_file, '-out', outprefix,]
        try:
            if self.verbose:
                print('build igblastp database: ', ' '.join(cmd))
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            if self.verbose:
                print("Command output:", result.stdout.strip())
            return outprefix
        except subprocess.CalledProcessError as e:
            if self.verbose:
                print(f"Command failed with error: {e}")
                print(f"Stderr: {e.stderr}")
        return None


    def igblastp(self, fa_file:str, db_path:str, outfile:str, outfmt:int=None):
        outfmt = 5 if outfmt is None else outfmt
        cmd = [
            'igblastp',
            '-query', fa_file,
            '-db', db_path,
            '-outfmt', str(outfmt),
            '-out', outfile,
        ]
        try:
            if self.verbose:
                print('blastp: ', ' '.join(cmd))
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            if self.verbose:
                print("Command output:", result.stdout.strip())
            return 'succeed'
        except subprocess.CalledProcessError as e:
            if self.verbose:
                print(f"Command failed with error: {e}")
                print(f"Stderr: {e.stderr}")
        return 'fail'
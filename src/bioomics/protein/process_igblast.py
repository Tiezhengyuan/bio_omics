'''
install blast: sudo apt install ncbi-blast+
'''
import os
import sys
import subprocess

from ..common.dir import Dir

class ProcessIgblast:

    def __init__(self, bin_dir:str, outdir:str=None, verbose:bool=False):
        # igblast bin
        self.bin_dir = bin_dir
        if self.bin_dir not in sys.path:
            sys.path.append(self.bin_dir)
        # outdir
        self.outdir = outdir
        if outdir:
            Dir(outdir).init_dir()
        # print information
        self.verbose = verbose

    def build_prot_db(self, fa_file:str, bin_dir:str=None):
        exe = 'makeblastdb' if bin_dir is None else \
            os.path.join(bin_dir, 'makeblastdb')
        outdir = self.outdir if self.outdir else os.path.dirname(fa_file)
        prefix = os.path.splitext(os.path.basename(fa_file))[0]
        outprefix = os.path.join(outdir, prefix)
        cmd = [exe, '-parse_seqids', '-dbtype', 'prot', \
            '-in', fa_file, '-out', outprefix,]
        try:
            if self.verbose:
                os.system(f"{exe} -version")
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
        # outfmt could be 3, 4, 7
        outfmt = 5 if outfmt is None else outfmt
        cmd = [
            os.path.join(self.bin_dir, 'igblastp'),
            '-germline_db_V', db_path,
            '-query', fa_file,
            '-outfmt', str(outfmt),
            '-out', outfile, ''
        ]
        try:
            if self.verbose:
                print('igblastp: ', ' '.join(cmd))
            result = subprocess.run(' '.join(cmd), shell=True, \
                capture_output=True, text=True, check=True)
            if self.verbose:
                print("Command output:", result.stdout.strip())
            return 'succeed'
        except subprocess.CalledProcessError as e:
            if self.verbose:
                print(f"Command failed with error: {e}")
                print(f"Stderr: {e.stderr}")
        return 'fail'
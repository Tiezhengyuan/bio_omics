'''
install blast: sudo apt install ncbi-blast+
'''
import os
import subprocess
from ..common import Dir

class ProcessBlast:

    def __init__(self, outdir:str=None, verbose:bool=False):
        self.outdir = outdir
        self.verbose = verbose

    def aln_iter(self, extension:str):
        for path in Dir(self.outdir).recursive_files():
            if path.endswith(extension):
                yield path

    def build_blastp_db(self, fa_file:str, bin_dir:str=None):
        db_type = 'prot'
        return self._blast_db(fa_file, db_type, bin_dir)

    def _blast_db(self, fa_file:str, db_type:str, bin_dir:str=None):
        outdir = self.outdir if os.path.isdir(self.outdir) else os.path.dirname(fa_file)
        Dir(outdir).init_dir()
        prefix = os.path.splitext(os.path.basename(fa_file))[0]
        cmd = [
            os.path.join(self.bin_dir, 'makeblastdb') if bin_dir else 'makeblastdb',
            '-in', fa_file,
            '-dbtype', db_type,
            '-out', os.path.join(outdir, prefix),
        ]
        try:
            if self.verbose:
                print('build blast database: ', ' '.join(cmd))
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            if self.verbose:
                print("Command output:", result.stdout.strip())
            return os.path.join(outdir, prefix)
        except subprocess.CalledProcessError as e:
            if self.verbose:
                print(f"Command failed with error: {e}")
                print(f"Stderr: {e.stderr}")
        return None


    def blastp(self, fa_file:str, db_path:str, outfile:str, outfmt:int=None):
        outfmt = 7 if outfmt is None else outfmt
        cmd = [
            'blastp',
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
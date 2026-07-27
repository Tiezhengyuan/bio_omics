'''
install blast
    wget https://ftp.debian.org/debian/pool/main/n/ncbi-blast+/ncbi-blast+_2.17.0+ds-7_amd64.deb
    sudo dpkg -i ncbi-blast+_2.17.0+ds-7_amd64.deb
'''
from pathlib import Path
import subprocess
from ..common import Dir

class ProcessBlast:

    def __init__(self, outdir:Path=None, verbose:bool=False):
        self.outdir = outdir
        self.verbose = verbose

    def aln_iter(self, extension:str):
        for path in Dir(self.outdir).recursive_files():
            if path.endswith(extension):
                yield path

    def _blast_db(self, fa_file:Path, db_type:str, bin_dir:Path=None):
        outdir = self.outdir if self.outdir.exists() else fa_file.parent
        outdir.mkdir(parents=True, exist_ok=True)
        db_path = outdir / fa_file.stem
        cmd = [
            str(bin_dir / 'makeblastdb') if bin_dir else 'makeblastdb',
            '-in', str(fa_file),
            '-dbtype', db_type,
            '-out', str(db_path),
        ]
        try:
            if self.verbose:
                print('build blast database: ', ' '.join(cmd))
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            if self.verbose:
                print("Command output:", result.stdout.strip())
            return db_path
        except subprocess.CalledProcessError as e:
            if self.verbose:
                print(f"Command failed with error: {e}")
                print(f"Stderr: {e.stderr}")
        return None

    def build_blastp_db(self, fa_file:Path, bin_dir:Path=None):
        db_type = 'prot'
        return self._blast_db(fa_file, db_type, bin_dir)

    def blastp(self, fa_file, db_path, outfile, outfmt:int=None):
        outfmt = 7 if outfmt is None else outfmt
        cmd = [
            'blastp',
            '-query', str(fa_file),
            '-db', str(db_path),
            '-outfmt', str(outfmt),
            '-out', str(outfile),
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
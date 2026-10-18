'''
install blast: sudo apt install ncbi-blast+
igblast: https://ftp.ncbi.nih.gov/blast/executables/igblast/release/
'''
from pathlib import Path
import sys
import subprocess


class ProcessIgblast:

    def __init__(self, igblast_dir:Path, verbose:bool=False):
        self.igblast_dir = Path(igblast_dir)
        self.verbose = verbose

    def build_pro_db(self, fa_file:Path, outdir:Path=None):
        fa_file = Path(fa_file)
        if outdir is None:
            outdir = fa_file.parent
        else:
            outdir.mkdir(parents=True, exist_ok=True)
        outprefix = str(outdir / fa_file.stem)
        
        cmd = [
            str(self.igblast_dir / 'bin' / 'makeblastdb'),
            '-parse_seqids',
            '-dbtype', 'prot',
            '-in', str(fa_file),
            '-out', outprefix,
        ]
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
        except FileNotFoundError as e:
            if self.verbose:
                print(f"Command not found: {e}")
        return None

    def igblastp(self, fa_file:Path, db_path:str, outfile:Path, outfmt:int=None):
        '''
        -outfmt <String> alignment view options:
            3 = Flat query-anchored, show identities,
            4 = Flat query-anchored, no identities,
            7 = Tabular with comment lines
            19 = Rearrangement summary report (AIRR format)
        '''
        outfmt = 5 if outfmt is None else outfmt
        cmd = [
            str(self.igblast_dir / 'bin' / 'igblastp'),
            '-query', str(fa_file),
            '-outfmt', str(outfmt),
            '-out', str(outfile),
            '-germline_db_V', db_path,
        ]
        try:
            if self.verbose:
                print('\n##igblastp: ', ' '.join(cmd))
            result = subprocess.run(
                cmd,
                cwd=str(self.igblast_dir),
                capture_output=True,
                text=True,
                check=True
            )
            if self.verbose:
                print("Command output:", result.stdout.strip())
            return 'succeed'
        except subprocess.CalledProcessError as e:
            if self.verbose:
                print(f"Command failed with error: {e}")
                print(f"Stderr: {e.stderr}")
        except FileNotFoundError as e:
            if self.verbose:
                print(f"Command not found: {e}")
        return 'fail'
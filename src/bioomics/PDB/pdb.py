'''
process PDB
1. renumber chains
'''
from Bio.PDB import PDBIO, PDBParser, MMCIFParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import gzip
from itertools import combinations
import json
from pathlib import Path

# Filter out the specific warning
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
warnings.simplefilter("ignore", PDBConstructionWarning)

from .atom_select import *
from .parse_structure import ParseStructure
from .pull_pdb import PullPdb
from ..protein.process_fasta import ProcessFasta


class Pdb(ParseStructure):

    def __init__(self, pdb_file:Path, outdir:Path=None, verbose:bool=False):
        super().__init__()
        # pdb_file, accept both str and Path
        self.pdb_file = Path(pdb_file)
        self.outdir = Path(outdir) if outdir else self.pdb_file.parent
        self.verbose = verbose
        # default objects should be updated
        self.annot = {}
        self.info = []
        # load structure
        self.load_structure()

    def load_pdb(self):
        self.parser = PDBParser(PERMISSIVE=1)
        self.structure = self.parser.get_structure(
            self.pdb_file.name, str(self.pdb_file)
        )

    def load_pdb_gz(self):
        self.parser = PDBParser(PERMISSIVE=1)
        with gzip.open(self.pdb_file, "rt") as f:
            self.structure = self.parser.get_structure(self.pdb_file.name, f)

    def load_cif(self):
        self.parser = MMCIFParser(QUIET=True, auth_chains=True)
        self.structure = self.parser.get_structure(
            self.pdb_file.name, str(self.pdb_file)
        )

    def load_cif_gz(self):
        self.parser = MMCIFParser(QUIET=True, auth_chains=False)
        with gzip.open(self.pdb_file, "rt") as f:
             self.structure = self.parser.get_structure(self.pdb_file.name, f)

    def load_structure(self):
        if not self.pdb_file.is_file():
            return None
        try:
            if self.pdb_file.name.endswith('.cif.gz'):
                self.load_cif_gz()
            elif self.pdb_file.name.endswith('.cif'):
                self.load_cif()
            elif self.pdb_file.name.endswith('.pdb.gz') or self.pdb_file.name.endswith('.ent.gz'):
                self.load_pdb_gz()
            else:
                self.load_pdb()

            if self.structure:
                self.structure_id = self.structure.header.get('idcode', '').upper()
                if self.verbose:
                    print(f"Successfully retrieve structure of {self.structure_id}")
        except Exception as e:
            print(f"Error: {self.pdb_file} is invalid. error={e}")
            self.structure = []
            self.structure_id = None

    def init_cascade_outdir(self, pdb_id:str=None):
        if pdb_id is None:
            pdb_id = self.structure_id
        if pdb_id:
            cascade_outdir = self.outdir / pdb_id[:2] / pdb_id
            cascade_outdir.mkdir(parents=True, exist_ok=True)
            if self.verbose:
                print('outputs dir: ', cascade_outdir)
            return cascade_outdir
        return None
    
    def init_chains_outdir(self, pdb_id:str=None):
        if pdb_id is None:
            pdb_id = self.structure_id
        if pdb_id:
            cascade_outdir = self.outdir / pdb_id[:2] / pdb_id / 'chains'
            cascade_outdir.mkdir(parents=True, exist_ok=True)
            if self.verbose:
                print('outputs dir: ', cascade_outdir)
            return cascade_outdir
        return None

    def header_lines(self):
        header_lines = []
        if self.pdb_file.suffix =='.gz':
            with gzip.open(self.pdb_file, 'rt') as f:
                for line in f:
                    if line.startswith(("HEADER", "TITLE", "COMPND", "SOURCE")):
                        header_lines.append(line)
        else:
            with open(self.pdb_file, 'r') as f:
                for line in f:
                    if line.startswith(("HEADER", "TITLE", "COMPND", "SOURCE")):
                        header_lines.append(line)
        return header_lines

    def save_annot(self):
        '''
        1. self.chains
        2. mapping of pdb - uniprot
        '''
        # resolution
        self.annot['resolution'] = self.structure.header.get('resolution')
            
        # integrate chains
        if self.chains:
            self.annot['chains'] = self.chains

        # integrate uniprot 
        uniprot = PullPdb(verbose=self.verbose).pdb_to_uniprot(self.structure_id)
        if uniprot:
            self.annot.update(uniprot)

        # save
        if self.annot:
            outfile = self.outdir / 'meta.json'
            if self.verbose:
                print(f"Save annotation of {self.structure_id}: {outfile}")
            with open(outfile, 'w') as f:
                json.dump(self.annot, f, indent=4)            

    def get_table(self):
        """
        get the score table from the bottom of a PDB
        return a list of lines
        """
        raw_table = []
        tag = False
        with open(self.pdb_file, 'r') as f:
            for line in f:
                if line.startswith("#BEGIN_POSE_ENERGIES_TABLE"):
                    tag = True
                if len(line) > 1 and tag is True:
                    raw_table.append(line)
        return raw_table

    def export_pdb(self, outfile_name:str, outdir:Path=None, has_table=False):
        outfile = outdir / outfile_name if outdir else self.pdb_file.parent / outfile_name
        if self.verbose:
            print(f"Try to create {outfile}...")
        with open(outfile, 'w') as f:
            # save atoms
            io=PDBIO()
            io.set_structure(self.structure)
            io.save(f)

            # save table
            if(has_table):
                raw_table = self.get_table()
                f.writelines(raw_table)
    
    def export_fasta(self, meta:defaultdict=None, overwrite:bool=False):
        '''
        retrieve sequence by chains from pdb
        '''
        if meta is None:
            meta = defaultdict(int)
        cascade_dir = self.init_cascade_outdir()
        outfile = cascade_dir / 'chain_sequences.faa'
        pdb_id = self.structure_id

        # prepare records
        unique_seq = []
        if overwrite or not outfile.is_file():
            records = []
            for item in self.parse_chains():
                chain_id, chain = item['chain_id'], item['chain']
                chain_seq = self._chain_seq(chain)
                # only export unique sequences to chain_sequences.faa
                if chain_seq and chain_seq not in unique_seq:
                    record = SeqRecord(
                        Seq(chain_seq),
                        id=chain_id,
                        description='',
                    )
                    records.append(record)
                    unique_seq.append(chain_seq)
            # save to fasta
            if records:
                status = ProcessFasta(outfile, self.verbose).to_fasta(records)
                meta[status] += 1
            else:
                meta['null'] += 1

        if outfile.is_file():
            return (pdb_id, outfile)
        return (pdb_id, None)

    # TODO: deprecated in the future
    def export_chain_faa(self, meta:defaultdict=None, overwrite:bool=False):
        '''
        retrieve sequence by chains from pdb
        '''
        db_records = []
        if meta is None:
            meta = defaultdict(int)
        pdb_id = self.structure_id
        cascade_dir = self.init_chains_outdir(pdb_id)
        # retrieve chains
        chains = self.parse_chains()
        for item in chains:
            chain_id, chain = item['chain_id'], item['chain']
            outfile = cascade_dir / f"{chain_id}.faa"
            if overwrite or not outfile.exists():
                chain_seq = self._chain_seq(chain)
                record = SeqRecord(
                    Seq(chain_seq),
                    id=chain_id,
                    description='',
                )
                if chain_seq:
                    status = ProcessFasta(outfile, self.verbose).to_fasta(record)
                    meta[status] += 1
                else:
                    meta['null'] += 1
            else:
                meta['skip'] += 1
            if outfile.is_file():
                rec = (chain_id, outfile)
                db_records.append(rec)
        return db_records

    # TODO: deprecated in the future

    def export_chain_pdb(self, meta:defaultdict=None, overwrite:bool=False):
        """
        export by single chain in .pdb
        Note: number chain pdb could be more than *.faa
            because some of structure are not protein
        """
        db_records = []
        if meta is None:
            meta = defaultdict(int)
        pdb_id = self.structure_id
        dir_endpoint = Path('data') / pdb_id[:2] / pdb_id / 'chain_pdb'
        cascade_dir = self.outdir / dir_endpoint
        cascade_dir.mkdir(parents=True, exist_ok=True)
        # retrieve chains
        chains = self.parse_chains()
        for item in chains:
            chain_id = item['chain_id']
            model_no = item['model_no']
            chain = item['chain']
            chain_seq = self._chain_seq(chain)
            # only consider limited models
            if len(chain_seq) > 0 and int(model_no) <= 10:
                file_name = f"{chain_id}.pdb"
                outfile = cascade_dir / file_name
                if overwrite or not outfile.is_file():
                    status = self.to_chain_pdb(model_no, [item['chain_no'],], outfile)
                    meta[status] += 1
                else:
                    meta['skip'] += 1
                relative_path = dir_endpoint / file_name
                rec = (chain_id, relative_path) if outfile.is_file() else (chain_id, None)
                db_records.append(rec)
        return db_records
    
    def export_combo_pdb(self, combo:int=None, meta:defaultdict=None, overwrite:bool=False):
        """
        export by two chains combinations
        """
        if combo is None:
            combo = 2
        if meta is None:
            meta = defaultdict(int)
        meta['pdb'] += 1

        chain_no_pool = self.get_chain_no()
        if not chain_no_pool:
            if self.verbose:
                print(f"Skip: No chains are detected in {self.pdb_file}.")
            meta['null'] += 1
            return None

        pdb_id = self.structure_id
        cascade_dir = self.init_chains_outdir(pdb_id)
        # export single/pairwise chain pdb
        outfiles, fail = [], []
        records = []
        for model_no, _pool in chain_no_pool.items():
            _pool.sort()
            for _ids in combinations(_pool, combo):
                chain_combo = '-'.join(list(_ids))
                outfile = cascade_dir / f"{pdb_id}_{model_no}_{chain_combo}.pdb"
                if overwrite is True or not outfile.is_file():
                    status = self.to_chain_pdb(model_no, list(_ids), outfile)
                    if status:
                        outfiles.append(outfile)
                    else:
                        fail.append(outfile)
                else:
                    meta['skip_build_faa'] += 1
                rec = (pdb_id, model_no, chain_combo, outfile)
                records.append(rec)
        if outfiles == []:
            rec = (pdb_id, None, None, None)
            records.append(rec)
        meta['succeessful_faa'] += len(outfiles)
        meta['fail_faa'] += len(fail)
        return records

    def to_chain_pdb(self, model_no:int, chain_no_pool:list, outfile):
        """
        Save chain while preserving original header information
        TODO: comfirm 1BAG A chain B is oligo
        """
        try:
            # Write output with header
            io = PDBIO(use_model_flag=model_no)
            io.set_structure(self.structure)
            # selector
            _ids = ''.join(chain_no_pool)
            atom_selector = ChainAtomSelect
            if '-' in _ids:
                atom_selector = LigandAtomSelect
            elif '+' in _ids:
                atom_selector = AllAtomSelect

            # Copy header from original file
            header_lines = self.header_lines()
            with open(outfile, "w") as f:
                f.writelines(header_lines)
                io.save(f, atom_selector(chain_no_pool))
            if self.verbose:
                print(f"Saved chain {chain_no_pool} to {outfile}")
            return True
        except Exception as e:
            print(f'Error: can not extract chain {chain_no_pool} from {self.structure_id}, error={e}')
        return False


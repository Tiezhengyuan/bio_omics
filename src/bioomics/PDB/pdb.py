'''
process PDB
1. renumber chains
'''
from collections import defaultdict
import gzip
from itertools import combinations
import json
import os
from Bio.PDB import PDBIO, PDBParser, MMCIFParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from biosequtils import Dir

# Filter out the specific warning
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
warnings.simplefilter("ignore", PDBConstructionWarning)

from .atom_select import *
from .parse_structure import ParseStructure
from ..constants.amino_acids import longer_names
from ..protein.process_fasta import ProcessFasta

# from .download import Download


class Pdb(ParseStructure):

    def __init__(self, pdb_file:str, outdir:str=None, verbose:bool=False):
        super().__init__()
        # pdb_file
        self.pdb_file = pdb_file
        self.outdir = outdir if outdir else os.path.dirname(pdb_file)
        self.verbose = verbose
        # default objects should be updated
        self.structure_id = None
        self.annot = {}
        self.info = []
        # load structure
        # self.load_structure()

    def get_parser(self):
        '''
        update self.pdb_file_name, self.parser
        '''
        self.pdb_file_name = os.path.basename(self.pdb_file)
        file_type = 'cif' if self.pdb_file_name.endswith('.cif') or \
            self.pdb_file_name.endswith('.cif.gzs')else 'pdb'

        self.parse = None
        if file_type == 'pdb':
            self.parser = PDBParser(PERMISSIVE=1)
        elif file_type == 'cif':
            self.parser = MMCIFParser(QUIET=True)

    def load_structure(self):
        """
        update self.structure, self.structure_id
        """
        # initialize parser
        self.get_parser()
        # load structure
        if self.pdb_file.endswith('.gz'):
            with gzip.open(self.pdb_file, 'rt') as gzf:
                self.structure = self.parser.get_structure(
                    self.pdb_file_name,
                    gzf
                )
        else:
            self.structure = self.parser.get_structure(
                self.pdb_file_name,
                self.pdb_file
            )
        if self.structure:
            self.structure_id = self.structure.header.get('idcode', '').upper()
            if self.verbose:
                print(f"Successfully retrieve structure of {self.structure_id}")

    def init_cascade_outdir(self, pdb_id:str=None):
        if pdb_id is None:
            pdb_id = self.structure_id
        if pdb_id:
            cascade_outdir = os.path.join(self.outdir, pdb_id[:2], pdb_id)
            Dir(cascade_outdir).init_dir()
            if self.verbose:
                print('outputs dir: ', cascade_outdir)
            return cascade_outdir
        return None
    
    def init_chains_outdir(self, pdb_id:str=None):
        if pdb_id is None:
            pdb_id = self.structure_id
        if pdb_id:
            cascade_outdir = os.path.join(self.outdir, pdb_id[:2], pdb_id, 'chains')
            Dir(cascade_outdir).init_dir()
            if self.verbose:
                print('outputs dir: ', cascade_outdir)
            return cascade_outdir
        return None

    def header_lines(self):
        header_lines = []
        if self.pdb_file.endswith('.gz'):
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
        uniprot = Download.pdb_to_uniprot(self.structure_id)
        if uniprot:
            self.annot.update(uniprot)

        # save
        if self.annot:
            outfile = os.path.join(self.outdir, 'meta.json')
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

    def export_pdb(self, outfile_name:str, outdir:str=None, has_table=False):
        indir = os.path.dirname(self.pdb_file)
        outfile = os.path.join(outdir, outfile_name) if outdir \
            else os.path.join(indir, outfile_name)
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
        cascade_dir = self.init_cascade_outdir()
        outfile = os.path.join(cascade_dir, 'chain_sequences.faa')
        pdb_id = self.structure_id

        # prepare records
        unique_seq = []
        if overwrite or (not os.path.isfile(outfile)):
            records = []
            for chain_id, _pdb_id, model_no, chain_no, chain in self.parse_chains():
                chain_seq = self._chain_seq(chain)
                # only export unique sequences to chain_sequences.faa
                if chain_seq not in unique_seq:
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

        if os.path.isfile(outfile):
            return (pdb_id, outfile)
        return None

      
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
        for chain_id, _pdb_id, model_no, chain_no, chain in chains:
            outfile = os.path.join(cascade_dir, f"{chain_id}.faa")
            if overwrite or (not os.path.isfile(outfile)):
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
            if os.path.isfile(outfile):
                db_records.append((chain_id, outfile))
        return db_records


    def export_chain_pdb(self, meta:defaultdict=None, overwrite:bool=False):
        """
        export by single chain, and two chains combinations
        """
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
            outprefix = os.path.join(cascade_dir, f"{pdb_id}_{model_no}_")
            _pool.sort()
            for combo in (1, 2,):
                for _ids in combinations(_pool, combo):
                    chain_combo = '-'.join(list(_ids))
                    outfile = outprefix + chain_combo + '.pdb'
                    if overwrite is True or not os.path.isfile(outfile):
                        status = self.export_chain_pdb(list(_ids), outfile)
                        if status:
                            outfiles.append(outfile)
                        else:
                            fail.append(outfile)
                    else:
                        meta['skip_build_faa'] += 1
                    rec = (pdb_id, model_no, chain_combo, outfile)
                    records.append(rec)
        meta['succeessful_faa'] += len(outfiles)
        meta['fail_faa'] += len(fail)
        return records


    def to_chain_pdb(self, chain_ids:list, outfile):
        """
        Save chain while preserving original header information
        """
        try:
            # Write output with header
            io = PDBIO()
            io.set_structure(self.structure)
            # selector
            _ids = ''.join(chain_ids)
            atom_selector = ChainAtomSelect
            if '-' in _ids:
                atom_selector = LigandAtomSelect
            elif '+' in _ids:
                atom_selector = AllAtomSelect

            # Copy header from original file
            header_lines = self.header_lines()
            with open(outfile, "w") as f:
                f.writelines(header_lines)
                io.save(f, atom_selector(chain_ids))
            if self.verbose:
                print(f"Saved chain {chain_ids} to {outfile}")
            return True
        except Exception as e:
            pass
        return False

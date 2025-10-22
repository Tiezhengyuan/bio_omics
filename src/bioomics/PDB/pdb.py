'''
process PDB
1. renumber chains
'''
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
    
    def retrieve_sequence(self) -> dict:
        records = {}
        for model in self.structure:
            model_no = model.id
            for chain in model:
                chain_no = chain.id.upper()
                seq = []
                for residue in chain:
                    if residue.resname in longer_names:
                        seq.append(longer_names[residue.resname])
                if seq:
                    key = (model_no, chain_no)
                    records[key] = SeqRecord(
                        Seq(''.join(seq)),
                        id=f"{self.structure_id}_{model_no}_{chain_no}",
                        description='',
                    )
        return records

    def export_fasta(self, outfile:str):
        '''
        retrieve sequence by chains from pdb
        '''
        # retrieve sequences
        sequences = self.retrieve_sequence()
        records = {}
        for _, record in sequences.items():
            seq = str(record.seq)
            if seq not in records:
                records[seq] = record
        if records:
            records = records.values()
            return ProcessFasta(outfile, self.verbose).to_fasta(records)
        return 'na'
       
    def export_chains_fasta(self, overwrite:bool=False):
        '''
        retrieve sequence by chains from pdb
        '''
        pdb_id = self.structure_id
        cascade_dir = self.init_chains_outdir(pdb_id)
        # retrieve sequences
        sequences = self.retrieve_sequence()
        seq_pool = []
        for key, record in sequences.items():
            model_no, chain_no = key
            seq = str(record.seq)
            if seq not in seq_pool:
                outfile = os.path.join(cascade_dir, f"{model_no}_{chain_no}.faa")
                if overwrite or not os.path.isfile(outfile):
                    ProcessFasta(outfile, self.verbose).to_fasta(record)
                seq_pool.append(seq)

    def export_chain_pdb(self, chain_ids:list, outdir:str=None, overwrite:bool=False):
        """
        Save chain while preserving original header information
        """
        file_name = '_'.join(chain_ids)
        if outdir is None or not os.path.isdir(outdir):
            outdir = self.outdir
        outfile = os.path.join(outdir, f"{self.structure_id}_{file_name}.pdb")
        if overwrite is False and os.path.isfile(outfile):
            if self.verbose:
                print(f"Skip: {outfile} exists")
            return None

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
        return outfile

    def split_by_chain(self, overwrite:bool=False):
        """
        export by single chain, and two chains combinations
        """
        chain_no_pool = self.get_chain_no()
        if not chain_no_pool:
            if self.verbose:
                print(f"Skip: No chains are detected in {self.pdb_file}.")
            return None
        chain_no_pool.sort()

        pdb_id = self.structure_id
        cascade_dir = self.init_chains_outdir(pdb_id)
        # export single/pairwise chain pdb
        outfiles = []
        for combo in (1, 2,):
            for _ids in combinations(chain_no_pool, combo):
                outfile = self.export_chain_pdb(list(_ids), \
                    cascade_dir, overwrite)
                outfiles.append(outfile)
        return outfiles
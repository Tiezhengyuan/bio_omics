'''
process PDB
1. renumber chains
'''
import gzip
import json
import os
import re
from pprint import pprint
from Bio.PDB import PDBIO, PDBParser, MMCIFParser, Select
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from biosequtils import Dir

# Filter out the specific warning
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
warnings.simplefilter("ignore", PDBConstructionWarning)


from .parse_structure import ParseStructure
from ..constants.amino_acids import longer_names
from ..common.process_pickle import ProcessPickle
# from .atom_select import ChainAtomSelect, LigandAtomSelect, AllAtomSelect
# from .download import Download

class NonHetSelect(Select):
    """Select only ATOM records (exclude HETATM)"""
    def accept_residue(self, residue):
        return residue.id[0] == " "  # Keep only standard residues

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
        # update self.outdir
        if pdb_id is None:
            pdb_id = self.structure_id
        if pdb_id:
            cascade_outdir = os.path.join(self.outdir, pdb_id[:2], pdb_id)
            Dir(cascade_outdir).init_dir()
            if self.verbose:
                print('outputs dir: ', cascade_outdir)
            return cascade_outdir
        return None

    def split_by_chain(self, chain_ids:list, file_name:str=None):
        """
            Save chain while preserving original header information
        """
        # Copy header from original file
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

        # save
        chain_ids = [re.sub(r'\-|\+', '', i) for i in chain_ids]
        file_name = '_'.join(chain_ids) if file_name is None else file_name
        outname = f"{self.structure_id}_{file_name}.pdb"
        outfile = os.path.join(self.outdir, outname)
        with open(outfile, "w") as f:
            f.writelines(header_lines)
            io.save(f, atom_selector(chain_ids))
        if self.verbose:
            print(f"Saved chain {chain_ids} to {outfile}")
        return outfile

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
    
    def retrieve_sequence(self):
        records = []
        for model in self.structure:
            for chain in model:
                seq = []
                for residue in chain:
                    if residue.resname in longer_names:
                        seq.append(longer_names[residue.resname])
                if seq:
                    rec = SeqRecord(
                        Seq(''.join(seq)),
                        id=f"{self.structure_id}_{model.id}_{chain.id}",
                        description='',
                    )
                    records.append(rec)
        return records

    def export_fasta(self, outfile):
        '''
        retrieve sequence by chains from pdb
        '''
        # load structure
        self.load_structure()
        # retrieve sequences
        records = self.retrieve_sequence()
        if records:
            try:
                with open(outfile, 'w') as f:
                    SeqIO.write(records, f, 'fasta')
                    if self.verbose:
                        print(f"Export chain sequences of {self.structure_id} to {outfile}.")
                return 'succeed'
            except Exception as e:
                if self.verbose:
                    print(e)
                return 'fail'
        return 'na'
        

    def build_pdb_obj(self):
        pdb_obj = Pdb(self.pdb_file)
        file_name = f"{pdb_obj.structure_id}_pdb.pickle"
        outfile = os.path.join(self.outdir, file_name)
        pp = ProcessPickle(outfile, self.verbose)
        return pp.save_pickle(pdb_obj)


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


from ..constants.amino_acids import longer_names
# from .atom_select import ChainAtomSelect, LigandAtomSelect, AllAtomSelect
# from .download import Download

class NonHetSelect(Select):
    """Select only ATOM records (exclude HETATM)"""
    def accept_residue(self, residue):
        return residue.id[0] == " "  # Keep only standard residues


class ParseStructure:
    def __init__(self):
        self.structure = None
        self.chains = []

    def get_header(self):
        header = {
            'resolution': self.structure.header['resolution'],
            'keywords': self.structure.header['keywords'],
            'compound': self.structure.header['compound'],
            'source': self.structure.header['source'],
        }
        return header

    def iter_res(self):
        '''
        iterate residues
        '''
        for model in self.structure:
            for chain in model:
                for res in chain:
                    yield model, chain, res

    def _chain_seq(self, chain):
        '''
        return AA seq given a chain object
        '''
        seq = []
        for residue in chain:
            if residue.resname in longer_names:
                seq.append(longer_names[residue.resname])
        return ''.join(seq)

    def get_chains_seq(self) -> dict:
        '''
        update self.chains
        '''
        res = {}
        for model in self.structure:
            for chain in model:
                chain_id = chain.id
                seq = self._chain_seq(chain)
                if seq:
                    res[chain_id] = seq
        return res

    def get_chains(self):
        '''
        update self.chains
        '''
        for model in self.structure:
            model_chains = {
                'model_id': model.id,
                'chains': [],
            }
            for chain in model:
                _chain = {
                    'chain_id': chain.id,
                    'AA residues': len(chain),
                    'seq': self._chain_seq(chain)
                }
                model_chains['chains'].append(_chain)
            self.chains.append(model_chains)
        return self.chains

    def parse_alpha_atoms(self, chain_id:str):
        '''
        Get all alpha atoms from Chain A
        '''
        for model in self.structure:
            for chain in model:
                if chain.id == chain_id:
                    ca_atoms = []
                    for residue in chain:
                        # Check if residue has a alpha atom
                        if "CA" in residue: 
                            ca_atoms.append(residue)
                    return ca_atoms
        return None

    def remove(self, chain_id:str, id_start:int, id_end:int):
        '''
        remove residues from a chain
        '''
        for model in self.structure:
            for chain in model:
                if chain.id == chain_id:
                    to_remove = [res for res in chain if id_start <= int(res.id[1]) <= id_end ]
                    for res in to_remove:
                        chain.detach_child(res.id)

    def trim(self, chain_id:str, id_start:int):
        '''
        remove residues from a chain
        '''
        for model in self.structure:
            for chain in model:
                if chain.id == chain_id:
                    to_remove = [res for res in chain if int(res.id[1]) >= id_start ]
                    for res in to_remove:
                        chain.detach_child(res.id)

    def renumber(self, start:int=None, norestart=True, preserve=False):
        '''
        renumber chains
        arguments:
            norestart: don't start renumbering at each chain
            preserve: preserve insertion code and heteroflags
        '''
        residue_id = 1 if start is None else start

        # process structure
        chain_id = ""
        for residue in self.structure.get_residues():
            chain = residue.get_parent()
            if chain_id != chain.get_id() and not norestart:
                chain_id = chain.get_id()
                residue_id = int(start)
            if preserve:
                hetero = residue.id[0]
                insert = residue.id[2]
                residue.id = (hetero, residue_id, insert)
            else:
                residue.id = (' ', residue_id, ' ')
            residue_id += 1
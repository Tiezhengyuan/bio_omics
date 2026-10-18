'''
create unit test cases for methods defined in the class SecondaryStructure
- input data of protein structures are in tests/data.
- test cases are put in the script test_secondary_structures.py
'''
from Bio import PDB
from scipy.spatial import cKDTree
from collections import defaultdict
from pathlib import Path

from .trajectory import Trajectory


class SecondaryStructure:
    def __init__(self, pdb_file:Path):
        self.pdb_file = str(pdb_file)
        self.load_pdb()
        
        # SS metrics
        self.valid_pairs = set() #clashes scores
        self.ss_counts = defaultdict(int)
        self.ss_interface_counts = defaultdict(int)
        self.plddts_interface = []
        self.plddts_ss = []
        
    def load_pdb(self):
        parser = PDB.PDBParser(QUIET=True)
        self.structure = parser.get_structure("protein", self.pdb_file)
        
    def calc_clash_score(self, threshold:float=2.4, only_ca:bool=True):
        '''clash score of chains'''
        # reset the state so repeated calls on the same instance do not
        # accumulate stale pairs from previous thresholds/atom selections
        self.valid_pairs = set()
        atoms = []
        atom_info = []
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        # Skip hydrogen atoms
                        if atom.element == "H":
                            continue
                        # If only calculating clashes based on CA atoms, skip non-CA atoms
                        if only_ca and atom.get_name() != "CA":
                            continue
                        atoms.append(atom.coord)
                        atom_info.append((chain.id, residue.id[1], atom.get_name(), atom.coord))

        # no atoms selected (empty structure): no clashes
        if not atoms:
            return 0
        tree = cKDTree(atoms)
        pairs = tree.query_pairs(threshold)

        for i, j in pairs:
            chain_i, res_i, name_i, coord_i = atom_info[i]
            chain_j, res_j, name_j, coord_j = atom_info[j]
            
            # Exclude clashes within the same residue
            if chain_i == chain_j and res_i == res_j:
                continue
            # Exclude directly sequential residues in the same chain for all atoms
            if chain_i == chain_j and abs(res_i - res_j) == 1:
                continue
            # If calculating sidechain clashes, only consider clashes between different chains
            if not only_ca and chain_i == chain_j:
                continue

            self.valid_pairs.add((i, j))

        return len(self.valid_pairs)

    def count_residues(self,
            dssp_path:Path,
            chain_id:str="B",
            atom_distance_cutoff:float=4.0,
        ) -> tuple:
        '''
        update self.ss_counts, self.ss_interface_counts,
            self.plddts_interface, and self.plddts_ss
        '''
        model = self.structure[0]
        # Calculate DSSP for the model
        dssp = PDB.DSSP(model, self.pdb_file, dssp=str(dssp_path))
       
        # Get chain and interacting residues once
        interacting_residues = Trajectory(self.pdb_file).hotspot_residues(chain_id, atom_distance_cutoff)
        interacting_residues = set(interacting_residues.keys())

        # Prepare to count residues
        chain = model[chain_id]
        for residue in chain:
            residue_id = residue.id[1]
            if (chain_id, residue_id) in dssp:
                ss = dssp[(chain_id, residue_id)][2]
                ss_type = "loop"
                if ss in ("H", "G", "I"):
                    ss_type = "helix"
                elif ss == "E":
                    ss_type = "sheet"

                self.ss_counts[ss_type] += 1

                if ss_type != "loop":
                    avg_plddt_ss = sum(atom.bfactor for atom in residue) / len(residue)
                    self.plddts_ss.append(avg_plddt_ss)

                if residue_id in interacting_residues:
                    self.ss_interface_counts[ss_type] += 1
                    avg_plddt_residue = sum(atom.bfactor for atom in residue) / len(residue)
                    self.plddts_interface.append(avg_plddt_residue)

    @staticmethod
    def calculate_percentages(total, helix, sheet):
        if total > 0:
            helix_percentage = round((helix / total) * 100, 2)
            sheet_percentage = round((sheet / total) * 100, 2)
            loop_percentage = round(((total - helix - sheet) / total) * 100, 2)
            return helix_percentage, sheet_percentage, loop_percentage
        return 0.0, 0.0, 0.0
    
    def total_percentages(self):
        '''percentage of total residues, helix, sheet'''
        total_residues = sum(self.ss_counts.values())        
        percentages = SecondaryStructure.calculate_percentages(
            total_residues,
            self.ss_counts["helix"],
            self.ss_counts["sheet"]
        )
        return percentages
        
    def interface_percentages(self):
        '''interface_percentage of total residues, helix, sheet'''
        total_interface_residues = sum(self.ss_interface_counts.values())
        interface_percentages = SecondaryStructure.calculate_percentages(
            total_interface_residues,
            self.ss_interface_counts["helix"],
            self.ss_interface_counts["sheet"]
        )
        return interface_percentages

    def i_plddt(self):
        if self.plddts_interface:
            i_plddt = round(sum(self.plddts_interface) / len(self.plddts_interface) / 100, 2)
            return i_plddt
        return 0
        
    def ss_plddt(self):
        if self.plddts_ss:
            ss_plddt = round(sum(self.plddts_ss) / len(self.plddts_ss) / 100, 2)
            return ss_plddt
        return 0

'''
prompt: 
    create unit test cases for methods defined in the class Trajectory
    - input data of protein structures are in tests/data.
    - test cases are put in the script test_trajactory.py
'''
from Bio import PDB
import json
import numpy as np
from pathlib import Path
import scipy

from bioomics import three_to_one_map


class Trajectory:
    parser = PDB.PDBParser(QUIET=True)

    def __init__(self, trajectory_pdb:Path):
        self.trajectory_pdb = trajectory_pdb
    
    def retrieve_traj(self, traj_chain:str='A'):
        # Parse the PDB files
        structure_trajectory = self.parser.get_structure("trajectory", self.trajectory_pdb)
        model_trajectory = structure_trajectory[0]
        # default: Extract chain or chain=A from trajectory_pdb
        chain_trajectory = model_trajectory[traj_chain]
        # Extract residues from chain A in trajectory_pdb
        residues_trajectory = [
            residue for residue in chain_trajectory 
            if PDB.Polypeptide.is_aa(residue, standard=True)
        ]
        return residues_trajectory

    def retriev_target(self, starting_pdb:Path, chain_ids_string):
        '''parse residues of target odb'''
        # Parse the PDB files
        structure_starting = self.parser.get_structure("starting", starting_pdb)
        model_starting = structure_starting[0]

        # Extract the specified chains from starting_pdb
        chain_ids = chain_ids_string.split(",")

        residues_starting = []
        for chain_id in chain_ids:
            chain_id = chain_id.strip()
            chain = model_starting[chain_id]
            for residue in chain:
                if PDB.Polypeptide.is_aa(residue, standard=True):
                    residues_starting.append(residue)
        return residues_starting
    
    @staticmethod
    def ca_residues(residues) -> list:
        '''Calpha atoms'''
        ca_residues =[]
        for residue in residues:
            if "CA" in residue:
                ca_residues.append(residue["CA"])
        return ca_residues
    
    def target_pdb_rmsd(self, starting_pdb:Path, chain_ids_string):
        '''RMSD between trajectory and target protein'''
        residues_starting = self.retriev_target(starting_pdb, chain_ids_string)
        residues_trajectory = self.retrieve_traj()

        # Ensure that both structures have the same number of residues
        min_length = min(len(residues_starting), len(residues_trajectory))
        residues_starting = residues_starting[:min_length]
        residues_trajectory = residues_trajectory[:min_length]

        # Collect Calpha atoms from the two sets of residues
        atoms_starting = Trajectory.ca_residues(residues_starting)
        atoms_trajectory = Trajectory.ca_residues(residues_trajectory)

        # Calculate RMSD using structural alignment: least-squares superpositions
        sup = PDB.Superimposer()
        sup.set_atoms(atoms_starting, atoms_trajectory)
        rmsd = sup.rms
        return round(rmsd, 2)

    def hotspot_residues(self, binder_chain="B", atom_distance_cutoff=4.0, traj_chain='A'):
        '''identify interacting residues at the binder interface'''
        # Parse the PDB file
        structure = self.parser.get_structure("complex", self.trajectory_pdb)

        # Get the specified chain
        binder_atoms = PDB.Selection.unfold_entities(structure[0][binder_chain], traj_chain)
        binder_coords = np.array([atom.coord for atom in binder_atoms])

        # Get atoms and coords for the target chain
        target_atoms = PDB.Selection.unfold_entities(structure[0][traj_chain], traj_chain)
        target_coords = np.array([atom.coord for atom in target_atoms])

        # Build KD trees for both chains
        binder_tree = scipy.spatial.cKDTree(binder_coords)
        target_tree = scipy.spatial.cKDTree(target_coords)

        # Prepare to collect interacting residues
        interacting_residues = {}

        # Query the tree for pairs of atoms within the distance cutoff
        pairs = binder_tree.query_ball_tree(target_tree, atom_distance_cutoff)

        # Process each binder atom's interactions
        for binder_idx, close_indices in enumerate(pairs):
            binder_residue = binder_atoms[binder_idx].get_parent()
            binder_resname = binder_residue.get_resname()

            # Convert three-letter code to single-letter code using the manual dictionary
            if binder_resname in three_to_one_map:
                aa_single_letter = three_to_one_map[binder_resname]
                for close_idx in close_indices:
                    target_residue = target_atoms[close_idx].get_parent()
                    interacting_residues[binder_residue.id[1]] = aa_single_letter

        return interacting_residues

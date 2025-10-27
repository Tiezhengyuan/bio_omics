'''
PAE predicted aligned errors
PAE measures the confidence in the relative position of two residues -
'''
import json
import numpy as np
import pandas as pd
from scipy import spatial
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt

from .pdb import Pdb

class CalculatePae:
    def __init__(self, pdb_file1, pdb_file2):
        # usually 1 is reference pdb, 2 is target pdb
        self.pdb_file1 = pdb_file1
        self.pdb_file2 = pdb_file2

    def calculate_pae(self, chain1=None, chain2=None):
        print(f"Try to calculatePAE: {self.pdb_file1} vs {self.pdb_file2}")
        
        # load coordinate of Ca atom
        coords1, residues1 = self.load_ca_coordinates(self.pdb_file1, chain1)
        coords2, residues2 = self.load_ca_coordinates(self.pdb_file2, chain2)
        # pairwise residues
        if len(coords1) != len(coords2):
            min_len = min(len(coords1), len(coords2))
            coords1 = coords1[:min_len]
            coords2 = coords2[:min_len]
            residues1 = residues1[:min_len]
            residues2 = residues2[:min_len]
            print(f"Empoly {min_len} residues to calculate PAE..")
        
        print("Try to align structure...")
        R = self.kabsch_rotation_matrix(coords1, coords2)
        # calculate center of the box
        ref_center = np.mean(coords1, axis=0)
        target_center = np.mean(coords2, axis=0)
        # apply rotation and move
        coords2_aligned = (coords2 - target_center) @ R.T + ref_center
        
        print("Try to calculate distance matrix...")
        # distance between each pair of observations from two different collections of inputs.
        dist_matrix1 = spatial.distance.cdist(coords1, coords1)
        dist_matrix2 = spatial.distance.cdist(coords2_aligned, coords2_aligned)
        
        print("Try to calculate PAE matrix...")
        pae_matrix = np.abs(dist_matrix1 - dist_matrix2)
        
        # pae matric should be symmetric
        pae_matrix = (pae_matrix + pae_matrix.T) / 2
        return pae_matrix
    

    def load_ca_coordinates(self, pdb_file, chain_id=None):
        pdb_obj = Pdb(pdb_file)
        pdb_obj.load_structure()

        coords = []
        residues_info = []
        for model, chain, residues in pdb_obj.ca_residue_iter():
            if chain_id is None or (chain_id and chain_id == chain.id):
                for residue in residues:
                    coords.append(residue['CA'].coord)
                    residues_info.append({
                        'chain': chain.id,
                        'residue_num': residue.id[1],
                        'residue_name': residue.resname,
                    })
        return np.array(coords), residues_info

    def kabsch_rotation_matrix(self, P, Q):
        """
        calculate optimal rotation matrix
        by minimizing RMSD between two pairs of points
        args: P reference coordinates (N,3)
        args: Q is target coordinates (N,3)
        return: R (3,3)
        """
        # center of coordinates
        P_centered = P - np.mean(P, axis=0)
        Q_centered = Q - np.mean(Q, axis=0)
        # covariance matrix
        H = P_centered.T @ Q_centered
        # SVD; singular value decomposition
        U, S, Vt = np.linalg.svd(H)
        # confirm right-hand coordinates
        d = np.sign(np.linalg.det(Vt.T @ U.T))
        D = np.eye(3)
        D[2, 2] = d
        # rotation matrix
        R = Vt.T @ D @ U.T
        return R

    def calculate_pae_statistics(self, pae_matrix):
        """
        PAE statistics
        """
        stats = {
            'matrix_shape': pae_matrix.shape,
            'global_mean': np.mean(pae_matrix),
            'global_median': np.median(pae_matrix),
            'global_std': np.std(pae_matrix),
            'global_min': np.min(pae_matrix),
            'global_max': np.max(pae_matrix),
            'residue_means': np.mean(pae_matrix, axis=1),
            'residue_medians': np.median(pae_matrix, axis=1),
            'confidence': {},
        }
        # confidence interval
        confidence = []
        for threshold in range(1, 51):
            # probability of PAE < threshold
            prob = np.sum(pae_matrix < threshold) / pae_matrix.size
            confidence.append(prob)
        stats['confidence'] = pd.DataFrame({
            'threshold': range(1, 51),
            'confidence': confidence,
        })
        return stats

    @staticmethod
    def heatmap(pae_matrix, figsize=None):
        if figsize is None:
            figsize = (6,5) 
        fig, ax = plt.subplots(figsize=figsize, layout='tight')
        colors = ['#006400', '#228B22', '#32CD32', '#90EE90', '#FFFFFF']
        cmap = mcolors.LinearSegmentedColormap.from_list('green_to_white', colors)
        # reverse y-axis
        im = ax.imshow(pae_matrix, cmap=cmap, origin='upper',
            vmin=0, vmax=min(30, np.max(pae_matrix)),
            extent=[1, pae_matrix.shape[1], 1, pae_matrix.shape[0]])
        plt.colorbar(im, ax=ax, label='PAE (Å)')
        ax.set_xlabel('Scored Residue')
        ax.set_ylabel('Aligned Residue')
        ax.set_title('Predicted Aligned Error (PAE)')
        return fig, ax

    @staticmethod
    def load_pae_json(infile):
        with open(infile, 'r') as f:
            pae = json.load(f)
        pae_matrix = np.array(pae[0]['predicted_aligned_error'])
        return pae_matrix
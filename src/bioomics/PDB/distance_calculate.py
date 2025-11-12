'''
distance between two atoms
'''
import numpy as np
import pandas as pd


import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array

from ..constants.amino_acids import EXTRA_LONGER_NAMES
from .analyze_complex import AnalyzeComplex

class DistanceCalculate:

    def __init__(self, pdb_file:str):
        # the pdb_file should contain exact two chains
        self.pdb_file = pdb_file
        # default
        self.dist1, self.dist2 = None, None
    
    def process(self, chain_combo:str):
        '''
        chain_combo: example 'A-C' stands for A and C
        '''
        chain1_id, chain2_id = chain_combo.split('-')
        dist1, dist2 = self.cal_dist(chain1_id, chain2_id)
        if dist1 is not None and dist2 is not None:
            if len(dist1) > 0 and len(dist2) > 0:
                motif1 = DistanceCalculate.parse_seq(dist1)
                motif2 = DistanceCalculate.parse_seq(dist2)
                return {
                    'chain1': chain1_id,
                    'dist1': dist1,
                    'motif1': motif1,
                    'chain2': chain2_id,
                    'dist2': dist2,
                    'motif2': motif2,
                }
    
    def cal_dist(self, chain1_id:str, chain2_id:str):
        try:
            # load MDA for simulation
            mdau = mda.Universe(self.pdb_file)

            # Select chains
            cond1 = "record_type hetatm" if chain1_id[-1] == '-' \
                else f"chainID {chain1_id}"
            cond2 = "record_type hetatm" if chain2_id[-1] == '-' \
                else f"chainID {chain2_id}"
            chain1 = mdau.select_atoms(cond1)
            chain2 = mdau.select_atoms(cond2)

            # Get the coordinates
            coords1 = chain1.positions
            coords2 = chain2.positions

            # Calculate the distances between each pair of atoms
            dists = distance_array(coords1, coords2)

            # parse atoms by the value of minimum distance
            pair_idx1 = np.apply_along_axis(np.argmin, 1, dists)
            df1 = pd.DataFrame({
                'res_no': chain1.resids,
                'res': chain1.resnames,
                'aa': [EXTRA_LONGER_NAMES.get(i, np.nan) for i in chain1.resnames],
                'atom': [a.name for a in chain1.atoms],
                'value': np.apply_along_axis(np.min, 1, dists),
                'pair_res_no': [chain2.resids[i] for i in pair_idx1],
                'pair_res': [chain2.resnames[i] for i in pair_idx1],
                'pair_atom': [chain2.atoms[i].name for i in pair_idx1],
            })
            df1['pair_aa'] = [EXTRA_LONGER_NAMES.get(i, np.nan) for i in df1['pair_res']]
            df1 = df1[df1['aa'].notna()]
            df1 = df1[df1['pair_aa'].notna()]
            self.dist1 = df1.groupby(['res_no']).agg('min').reset_index()
            
            # 
            pair_idx2 = np.apply_along_axis(np.argmin, 0, dists)
            df2 = pd.DataFrame({
                'res_no': chain2.resids,
                'res': chain2.resnames,
                'aa': [EXTRA_LONGER_NAMES.get(i, np.nan) for i in chain2.resnames],
                'atom': [a.name for a in chain2.atoms],
                'value': np.apply_along_axis(np.min, 0, dists),
                'pair_res_no': [chain1.resids[i] for i in pair_idx2],
                'pair_res': [chain1.resnames[i] for i in pair_idx2],
                'pair_atom': [chain1.atoms[i].name for i in pair_idx2],
            })
            df2['pair_aa'] = [EXTRA_LONGER_NAMES.get(i, np.nan) for i in df2['pair_res']]
            df2 = df2[df2['aa'].notna()]
            df2 = df2[df2['pair_aa'].notna()]
            self.dist2 = df2.groupby(['res_no']).agg('min').reset_index()
            return self.dist1, self.dist2
        except Exception as e:
            print(f"ERROR: calculate distance, pdb={self.pdb_file}, err={e}")
        return None, None

    @staticmethod
    def parse_seq(dist, max_span:int=None, cutoff:float=None):
        '''
        present contacts interface. In default:
            - max span of AA is 10
            - max distance between atoms is 5 A
        '''
        if not isinstance(dist, pd.DataFrame):
            return None
        max_span = 10 if max_span is None else int(max_span)
        cutoff = 5 if cutoff is None else float(cutoff)        

        # filter dist: values less than 5 in default
        sig = dist[dist['value']<=cutoff]
        # group significant point
        group = AnalyzeComplex.group_significant_points(sig, max_span)

        # convert to seq
        motif = []
        for g in group:
            start, end = g[0][0], g[-1][0]
            sub = dist.loc[start:end,:]
            # seq should be continuous seq, but pair_aa may not be continuous
            motif.append({
                'seq': ''.join(sub['aa']),
                'pair_aa': ''.join(sub['pair_aa']),
                'sig_res': len(g),
                'start': g[0][0],
                'end': g[-1][0],
                'cutoff': cutoff,
                'max_span': max_span,
            })
        return motif

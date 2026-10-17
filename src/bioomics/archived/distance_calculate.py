'''
distance between two atoms
'''
import itertools
from pathlib import Path
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array

from ..constants.amino_acids import EXTRA_LONGER_NAMES
from .distance_analyze import DistanceAnalyze
from .pdb import Pdb

class DistanceCalculate:

    def __init__(self, pdb_file:Path, args:dict):
        '''
        pdb_file: pdb should contain exact two chains
        dist_cutoff: 10A
        max_span: max 10 amino acids between two contacts
        '''
        self.pdb_file = pdb_file
        self.verbose = args.get('verbose', False)
        self.dist_cutoff = float(args.get('dist_cutoff', 5))
        self.max_span = int(args.get('max_span', 10))
        self.dist1, self.dist2 = None, None
    
    def cal_dist(self, chain1_id:str, chain2_id:str):
        try:
            # load MDA for simulation
            mdau = mda.Universe(str(self.pdb_file))

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
    
    def process(self, chain_combo:str) -> dict:
        '''
        chain_combo: example 'A-C' stands for A and C
        '''
        chain1_id, chain2_id = chain_combo.split('-')
        dist1, dist2 = self.cal_dist(chain1_id, chain2_id)
        if dist1 is not None and dist2 is not None:
            if len(dist1) > 0 and len(dist2) > 0:
                motif1 = DistanceAnalyze.parse_seq(dist1, self.dist_cutoff, self.max_span)
                motif2 = DistanceAnalyze.parse_seq(dist2, self.dist_cutoff, self.max_span)
                if self.verbose:
                    n1, n2 = len(motif1), len(motif2)
                    print(f"input: {self.pdb_file}, contacts: {n1} vs {n2}")
                return {
                    'chain_combo': chain_combo,
                    'chain1': chain1_id,
                    'dist1': dist1,
                    'motif1': motif1,
                    'chain2': chain2_id,
                    'dist2': dist2,
                    'motif2': motif2,
                }
        return {}
    
    def __call__(self) -> list:
        res = []

        p = Pdb(self.pdb_file)
        p.load_structure()
        chain_no_pool = p.get_model_chain_no(model_no=0)
        if len(chain_no_pool) > 1:
            for (chain1_id, chain2_id) in itertools.combinations(chain_no_pool, 2):
                dist1, dist2 = self.cal_dist(chain1_id, chain2_id)
                if dist1 is not None and dist2 is not None:
                    if len(dist1) > 0 and len(dist2) > 0:
                        motif1 = DistanceAnalyze.parse_seq(dist1, self.dist_cutoff, self.max_span)
                        motif2 = DistanceAnalyze.parse_seq(dist2, self.dist_cutoff, self.max_span)
                        res.append({
                            'dist_cutoff': self.dist_cutoff,
                            'max_span': self.max_span,
                            'chain_combo': f"{chain1_id}-{chain2_id}",
                            'chain1': chain1_id,
                            'dist1': dist1,
                            'motif1': motif1,
                            'chain2': chain2_id,
                            'dist2': dist2,
                            'motif2': motif2,
                        })
        return res
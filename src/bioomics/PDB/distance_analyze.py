'''
'''
import os
import numpy as np
import pandas as pd

from ..common import ProcessPickle
from .analyze_complex import AnalyzeComplex

class DistanceAnalyze:
    def __init__(self, args):
        '''
        present contacts interface. In default:
            - max span of AA >=1, default 10
            - max distance between atoms is 5 A
        '''
        self.pfile = args.get('pickle_file', None)
        self.data = ProcessPickle(self.pfile).load_pickle({})
        self.max_span = int(args.get('max_span', 10))
        self.dist_cutoff = float(args.get('dist_cutoff', 5))
        self.verbose = args.get('verbose', False)

    def get_distance(self, pair:bool=False, rank_idx:int=0):
        col = 'dist2' if pair is True else 'dist1'
        dist = self.data[col].sort_values('value').iloc[rank_idx]
        dist_value = round(dist['value'], 2)
        res1 = f"{dist['aa']}-{dist['res_no']}"
        res2 = f"{dist['pair_aa']}-{dist['pair_res_no']}"
        return dist_value, res1, res2

    def distance_threshold(self):
        dist1 = np.round(self.data['dist1']['value'].sort_values())
        dist2 = np.round(self.data['dist2']['value'].sort_values())
        values = []
        for i in [0, 4, 9]:
            val1 = dist1[i] if i < len(dist1) else None
            values.append(val1)
            val2 = dist2[i] if i < len(dist2) else None
            values.append(val2)
        return values

    def retrieve_motif(self, min_seq_len:int=1):
        '''
        the motif is binding region
        '''
        motifs = []
        try:
            motif1 = self.parse_seq(self.data.get('dist1'))
            for m in motif1:
                if len(m['seq']) >= min_seq_len:
                    m['chain_no'] = self.data['chain1']
                    m['pair_chain_no'] = self.data['chain2']
                    m['pdb_id'] = self.data['pdb_id']
                    m['combo_id'] = self.data['combo_id']
                    motifs.append(m)
            motif2 = self.parse_seq(self.data.get('dist2'))
            for m in motif2:
                if len(m['seq']) >= min_seq_len:
                    m['chain_no'] = self.data['chain2']
                    m['pair_chain_no'] = self.data['chain1']
                    m['pdb_id'] = self.data['pdb_id']
                    m['combo_id'] = self.data['combo_id']
                    motifs.append(m)
        except Exception as e:
            print(self.pfile, e)
        return motifs

    def parse_seq(self, dist):
        if not isinstance(dist, pd.DataFrame):
            return None
        # filter dist: values less than 5 in default
        sig = dist[dist['value']<=self.dist_cutoff]
        # group significant point
        group = AnalyzeComplex.group_significant_points(sig, self.max_span)

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
                'dist_cutoff': self.dist_cutoff,
                'max_span': self.max_span,
            })
        return motif

    def get_chain_data(self, type_name:str, chain_len:int=None):
            df1, df2 = None, None
            chain_id = (self.data['chains'][0][0], self.data['chains'][1][0])
            if self.data['type'] == type_name:
                df1 = self.data['dist']['dist1']
                df2 = self.data['dist']['dist2']
            if (df1 is not None) and (df2 is not None):
                if chain_len is None or (chain_len and chain_len == len(df1)):
                    return chain_id, df1, df2
    
    def get_receptor_data(self, chain_len:int=None):
        chain_id = (self.data['chains'][0][0], self.data['chains'][1][0])
        df = self.data['dist']['dist2']
        if chain_len is None or (chain_len and chain_len == len(df)):
            return chain_id, df

    def print_min_dist(self):
        g1, g2 = self.data['dist1'], self.data['dist2']
        print('---'*10)
        sig1 = g1[g1['value']<=self.dist_cutoff]
        print(sig1)
        print('---'*10)
        sig2 = g2[g2['value']<=self.dist_cutoff]
        print(sig2)



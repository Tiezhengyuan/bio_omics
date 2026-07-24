'''
'''
import numpy as np
import pandas as pd

from .analyze_complex import AnalyzeComplex

class DistanceAnalyze:

    @staticmethod
    def get_distance(data, pair:bool=False, rank_idx:int=0):
        col = 'dist2' if pair is True else 'dist1'
        dist = data[col].sort_values('value').iloc[rank_idx]
        dist_value = round(dist['value'], 2)
        res1 = f"{dist['aa']}-{dist['res_no']}"
        res2 = f"{dist['pair_aa']}-{dist['pair_res_no']}"
        return dist_value, res1, res2

    @staticmethod
    def distance_threshold(data:dict):
        dist1 = np.round(data['dist1']['value'].sort_values())
        dist2 = np.round(data['dist2']['value'].sort_values())
        values = []
        for i in [0, 4, 9]:
            val1 = dist1[i] if i < len(dist1) else None
            values.append(val1)
            val2 = dist2[i] if i < len(dist2) else None
            values.append(val2)
        return values

    @staticmethod
    def parse_seq(dist, dist_cutoff:float=5, max_span:int=10) -> list:
        if not isinstance(dist, pd.DataFrame):
            return None
        # filter dist: values less than 5 in default
        sig = dist[dist['value']<=dist_cutoff]
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
                'dist_cutoff': dist_cutoff,
                'max_span': max_span,
            })
        return motif

    @staticmethod
    def retrieve_motif(data, min_seq_len:int=1):
        '''
        the motif is binding region
        '''
        motifs = []
        try:
            motif1 = DistanceAnalyze.parse_seq(data.get('dist1'))
            for m in motif1:
                if len(m['seq']) >= min_seq_len:
                    m['chain_no'] = data['chain1']
                    m['pair_chain_no'] = data['chain2']
                    m['pdb_id'] = data['pdb_id']
                    m['combo_id'] = data['combo_id']
                    motifs.append(m)
            motif2 = DistanceAnalyze.parse_seq(data.get('dist2'))
            for m in motif2:
                if len(m['seq']) >= min_seq_len:
                    m['chain_no'] = data['chain2']
                    m['pair_chain_no'] = data['chain1']
                    m['pdb_id'] = data['pdb_id']
                    m['combo_id'] = data['combo_id']
                    motifs.append(m)
        except Exception as e:
            print('Error:', e)
        return motifs
    
    @staticmethod
    def get_chain_data(data, type_name:str, chain_len:int=None):
            df1, df2 = None, None
            chain_id = (data['chains'][0][0], data['chains'][1][0])
            if data['type'] == type_name:
                df1 = data['dist']['dist1']
                df2 = data['dist']['dist2']
            if (df1 is not None) and (df2 is not None):
                if chain_len is None or (chain_len and chain_len == len(df1)):
                    return chain_id, df1, df2
    
    @staticmethod
    def get_receptor_data(data, chain_len:int=None):
        chain_id = (data['chains'][0][0], data['chains'][1][0])
        df = data['dist']['dist2']
        if chain_len is None or (chain_len and chain_len == len(df)):
            return chain_id, df

    @staticmethod
    def print_min_dist(data:dict, dist_cutoff:float=5):
        g1, g2 = data['dist1'], data['dist2']
        print('---'*10)
        sig1 = g1[g1['value']<=dist_cutoff]
        print(sig1)
        print('---'*10)
        sig2 = g2[g2['value']<=dist_cutoff]
        print(sig2)



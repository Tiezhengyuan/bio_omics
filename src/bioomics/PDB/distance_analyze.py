'''
'''
import os
import re


# from .analyze_complex import AnalyzeComplex
from .distance_calculate import DistanceCalculate
from ..common import ProcessPickle

class DistanceAnalyze:
    def __init__(self, data, max_span:int=None, cutoff:float=None, verbose:bool=None):
        self.data = data
        self.max_span = int(max_span) if max_span else 30
        self.verbose = verbose if verbose else False
        self.cutoff = float(cutoff) if cutoff else 5

    def iter_dist(self):
        '''
        return distance df
        _data is dict, keys are chain1, dist1, motif1
            chain2, dist2, motif2, and chain_combo
        '''
        for _data in self.data:
            dist1 = _data['dist1']
            dist2 = _data['dist2']
            yield dist1, dist2


    def print_min_dist(self, ix:int):
        '''
        '''
        dist_data = self.data[ix]['dist']
        # print(dist_data['type'])

        g1, g2 = dist_data['dist1'], dist_data['dist2']
        print('---'*10)
        sig1 = g1[g1['value']<=self.cutoff]
        print(sig1)
        print('---'*10)
        sig2 = g2[g2['value']<=self.cutoff]
        print(sig2)

    def retrieve_motif(self, outdir:str=None):
        '''
        the motif is binding region
        '''
        pairs = []
        for _data in self.data:
            motif1 = DistanceCalculate.parse_seq(_data['dist1'], \
                self.max_span, self.cutoff)
            motif2 = DistanceCalculate.parse_seq(_data['dist2'], \
                self.max_span, self.cutoff)
            if motif1 and motif2:
                pairs.append((_data['chain1'], motif1, _data['chain2'], motif2))
        # export
        if os.path.isdir(outdir):
            cond = f"dist-{self.cutoff}_span-{self.max_span}"
            outfile = os.path.join(self.outdir, f'seq_dist_{cond}.p')
            if self.verbose:
                print(f"Sequence-Distance data: {outfile}.")
            ProcessPickle(outfile).save_pickle(pairs)
            return pairs, outfile
        return pairs, None


    def get_chain_data(self, type_name:str, chain_len:int=None):
        res = {}
        for _data in self.data:
            df1, df2 = None, None
            chain_id = (_data['chains'][0][0], _data['chains'][1][0])
            if _data['type'] == type_name:
                df1 = _data['dist']['dist1']
                df2 = _data['dist']['dist2']
            if (df1 is not None) and (df2 is not None):
                if chain_len is None or (chain_len and chain_len == len(df1)):
                    res[chain_id] = (df1, df2)
        return res
    
    def get_receptor_data(self, chain_len:int=None):
        # 
        res = {}
        for _data in self.data:
            chain_id = (_data['chains'][0][0], _data['chains'][1][0])
            df = _data['dist']['dist2']
            if chain_len is None or (chain_len and chain_len == len(df)):
                res[chain_id] = df
        return res
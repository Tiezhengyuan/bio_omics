'''
process distance data
'''
import pandas as pd
import os

from .distance_analyze import DistanceAnalyze

class DistanceProcess:
    def __init__(self, data_dir:str, span:int, cutoff:float):
        self.data_dir = data_dir
        self.span = span
        self.cutoff = cutoff

    def retrieve_interface_seq(self, dist_meta, min_seq_len:int=6):
        res = []
        for row in dist_meta:
            args = {
                'pickle_file': os.path.join(self.data_dir, row['relative_pkl']),
                'max_span': self.span,
                'dist_cutoff': self.cutoff,
            }
            da = DistanceAnalyze(args)
            _motifs = da.retrieve_motif(min_seq_len=min_seq_len)
            res.extend(_motifs)
        df = pd.DataFrame(res)
        print('pdb', df['pdb_id'].nunique())
        print('combo', df['combo_id'].nunique())
        return df
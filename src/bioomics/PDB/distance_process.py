'''
process distance data
'''
import pandas as pd
import os

from .distance_analyze import DistanceAnalyze

class DistanceProcess:
    def __init__(self, args:dict):
        self.data_dir = args['data_dir']
        self.span = args.get('max_span', 10)
        self.cutoff = args.get('dist_cutoff', 5)
        self.min_seq_len = args.get('min_seq_len', 6)
        self.verbose = args.get('verbose', False)

    def retrieve_interface_seq(self, dist_meta):
        res = []
        print('Number of distance', len(dist_meta))
        for i, row in enumerate(dist_meta):
            pfile = os.path.join(self.data_dir, row['relative_pkl'])
            args = {
                'pickle_file': pfile,
                'max_span': self.span,
                'dist_cutoff': self.cutoff,
            }
            da = DistanceAnalyze(args)
            _motifs = da.retrieve_motif(self.min_seq_len)
            res.extend(_motifs)
        df = pd.DataFrame(res)
        return df
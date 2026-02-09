'''
process distance data
'''
import pandas as pd
import os

from .distance_analyze import DistanceAnalyze
from ..common import ProcessPickle

class DistanceProcess:
    def __init__(self, args:dict):
        self.data_dir = args['data_dir']
        self.span = args.get('max_span', 10)
        self.cutoff = args.get('dist_cutoff', 5)
        self.min_seq_len = args.get('min_seq_len', 6)
        self.verbose = args.get('verbose', False)

    def count_contacts(self, dist_meta, cutoff_pool):
        print('Number of distance', len(dist_meta))
        res = {i:[] for i in cutoff_pool}
        for row in dist_meta:
            pfile = os.path.join(self.data_dir, row['relative_pkl'])
            try:
                data = ProcessPickle(pfile).load_pickle()
                dist1 = data['dist1']
                dist2 = data['dist2']
                for cutoff in cutoff_pool:
                    num_ca1 =  len(dist1[dist1['value'] < cutoff])
                    num_ca2 = len(dist2[dist2['value'] < cutoff])
                    if num_ca1 or num_ca2:
                        rec = {
                            'pdb_id': row['pdb_id'],
                            'combo_id': row['combo_id'],
                            'chain_combo': row['chain_combo'],
                            'num_ca1': num_ca1,
                            'num_ca2': num_ca2,
                            'num_ca': (num_ca1 + num_ca2)/2,
                        }
                        res[cutoff].append(rec)
            except Exception as e:
                pass
        return res

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


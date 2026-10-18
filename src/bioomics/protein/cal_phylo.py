'''
calculate phylogentic distance
MSA is fasta format, determined by MUSCLE
'''
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
import pandas as pd

from ..common import Dir

class CalPhylo:

    def __init__(self, indir:str=None):
        self.indir = indir
        self.methods = ['identity', 'blosum62']
        self.data = {}

    def calc_scan(self):
        res = {i:[] for i in self.methods}
        file_iter = Dir(self.indir).recursive_files()
        for path in file_iter:
            if path.endswith('.afa'):
                scores = self.msa_distance(path)
                for method, dist in scores.items():
                    res[method].append(dist)
        # to df
        for method in self.methods:
            res[method] = pd.concat(res[method], axis=1).T \
                if res[method] else pd.DataFrame()
        return res

    def msa_distance(self, msa_fa_file:str):
        '''
        distance of each sequence to the first one in the MSA
        return: {method: pd.Series(index=names[1:], values=distances)}
        '''
        self.data['aln_file'] = msa_fa_file
        # Load multiple sequence alignment
        alignment = AlignIO.read(msa_fa_file, "fasta")
        # Calculate distance matrix
        scores = {}
        for method in self.methods:
            try:
                calculator = DistanceCalculator(method)
                dm = calculator.get_distance(alignment)
                dist = pd.DataFrame(dm.matrix)[0][1:]
                dist.name = dm.names[0]
                dist.index = dm.names[1:]
                scores[method] = dist
            except Exception as e:
                pass
        if scores:
            self.data['distance'] = scores
        return scores

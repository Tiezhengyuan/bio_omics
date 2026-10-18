import re
import numpy as np
import pandas as pd
from collections import Counter

from ..constants.constants import *

class EncodeAA:
    def __init__(self):
        self.hydro_ph7 = {s:PROPERTY[s]['hydrophobicity_ph7'] for s in PROPERTY}
        self.hydro = {s:PROPERTY[s]['hydrophobicity'] for s in PROPERTY}
        self.polar = {s:PROPERTY[s]['polarity'] for s in PROPERTY}
        self.pz = {s:PROPERTY[s]['polarizability'] for s in PROPERTY}
        self.vdw = {s:PROPERTY[s]['van_der_Waals_volume'] for s in PROPERTY}
        self.aa = list(PROPERTY)
        self.aa2 = self.aa_2()
        # self.aa3 = self.aa_3()
        self.smiles = {a:' '.join(list(b)) for a,b in SMILES_1.items()}
        # for a,b in SMILES_1.items():
        #     b = re.sub(r'\(=O\)|\(=N\)|\(=CN2\)|\(|\)', '', b)
        #     self.smiles[a] = re.sub(r'\d|\=', '', b)

    def aa_2(self):
        '''
        2 AA
        '''
        res = []
        for a in self.aa:
            for b in self.aa:
                res.append(a+b)
        return res

    def aa_3(self):
        '''
        3 AA
        '''
        res = []
        for a in self.aa:
            for b in self.aa:
                for c in self.aa:
                    res.append(a+b+c)
        return res

    def aa_smiles(self, seq:str):
        '''
        represent aa by smiles
        '''
        # only consider 20 AA, remove X unknown residues
        if not set(seq).difference(set(self.smiles)):
            seq = seq.replace('G', '')
            smiles = [self.smiles[s] for s in seq if s in self.smiles]
            _s = str(' '.join(smiles))
            _s = re.sub(r'^[A-Z| ]', '', _s)
            return _s

    def physcial_chemical(self, seq:str):
        res = [
            ['sequencing is',] + list(seq),
            ['hydrophobicity at PH7 is',] + [str(self.hydro_ph7.get(s, 0)) for s in seq],
            ['hydrophobicity is',] + [str(self.hydro.get(s, 0)) for s in seq],
            ['polarity is',] + [str(self.polar.get(s, 0)) for s in seq],
            ['polarizability is',] + [str(self.pz.get(s, 0)) for s in seq],
            ['van der Waals_volume is',] + [str(self.vdw.get(s, 0)) for s in seq],
        ]
        res = [' '.join(i) for i in res]
        if not res:
            print(seq)
        return res

    def hydrophobicity_ph7(self, seq:str):
        res = [self.hydro_ph7.get(s, 0) for s in seq]
        return res

    def hydrophobicity(self, seq:str):
        res = [self.hydro.get(s, 0) for s in seq]
        return res

    def polarity(self, seq:str):
        res = [self.polar.get(s, 0) for s in seq]
        return res

    def polarizability(self, seq:str):
        res = [self.pz.get(s, 0) for s in seq]
        return res

    def van_der_Waals_volume(self, seq:str):
        res = [self.vdw.get(s, 0) for s in seq]
        return res

    def freq_1aa(self, seq:str) -> dict:
        '''
        frequency of 1-AA
        '''
        counts = Counter(seq)
        n = len(seq)
        freq = {k:round(v/n,3) for k,v in counts.items()}
        return freq
    
    def freq_2aa(self, seq:str) -> dict:
        '''
        frequency of 2-AA
        '''
        n = len(seq)
        freq = {}
        for i in range(0, len(seq)-1):
            aa = seq[i:i+2]
            key = f"freq_{aa}"
            if key not in freq:
                freq[key] = 1
            else:
                freq[key] += 1
        freq = {k:round(v/n, 3) for k,v in freq.items()}
        return freq

#################
    def physical_chemical_text(self, seq:str) -> str:
        '''
        physical chemical properties given a epitope sequence
        '''
        # physical-chemical properties
        # names = ['hydrophobicity at PH7', 
        # 'hydrophobicity', 'polarity', \
        #     'polarizability', 'van der Waals volume']
        phyche = [
            # self.hydrophobicity_ph7(seq),
            # self.hydrophobicity(seq),
            # self.polarity(seq),
            # self.polarizability(seq),
            self.van_der_Waals_volume(seq),
        ]
        res = [' '.join(list(seq)), ] + [str(int(np.mean(i))) for i in phyche]
        res = ' | '.join(res)
        return res

    def physical_chemical(self, seq:str, label:int=None) -> pd.Series:
        '''
        physical chemical properties given a epitope sequence
        '''
        # physical-chemical properties
        names = ['hydro_ph7', 'hydro', 'polar', 'polar_stab', 'VDWV']
        phyche = [
            self.hydrophobicity_ph7(seq),
            self.hydrophobicity(seq),
            self.polarity(seq),
            self.polarizability(seq),
            self.van_der_Waals_volume(seq),
        ]
        phyche_mean = pd.Series([np.mean(i) for i in phyche],
            index=[f"mean_{i}" for i in names])
        phyche_median = pd.Series([np.median(i) for i in phyche],
            index=[f"median_{i}" for i in names])
        phyche_var = pd.Series([np.var(i) for i in phyche],
            index=[f"variance_{i}" for i in names])

        # outcome
        res = pd.Series([seq, label], index=['seq', 'label'])
        # round numeric stats before merging with string values
        stats = pd.concat([phyche_mean, phyche_median, phyche_var])
        stats = stats.round(decimals=3)
        res = pd.concat([res, stats])
        return res

    def frequency_aa(self, seq:str, label:int=None) -> list:
        '''
        frequency of amino acids
        '''
        res = [seq, label, ]
        # single
        freq1a = self.freq_1aa(seq)
        for k in self.aa:
            v = freq1a.get(k, 0)
            res.append(str(v))
        # two AA
        freq2a = self.freq_2aa(seq)
        for k in self.aa2:
            v = freq2a.get(f'freq_{k}', 0)
            res.append(str(v))
        return res

    def existing(self, seq:str, label:int=None) -> list:
        '''
        if an AA is existing
        binary: 1 or 0
        '''
        res = [seq, label]
        for pool in [self.aa, self.aa2]:
            for i in pool:
                if i in seq:
                    res.append(str(1))
                else:
                    res.append(str(0))
        return res

    def vector_1d(self, seq:str, label:int=None) -> np.array:
        res = self.hydrophobicity_ph7(seq) + self.hydrophobicity(seq) + \
                self.polarity(seq) + self.polarizability(seq) + \
                self.van_der_Waals_volume(seq)
        res = np.array(res)
        if label is not None:
            res = np.append(res, label)
        return res
    
    def vector_2d(self, seq:str) -> np.array:
        # one row per residue, one column per property
        res = self.vector_1d(seq)
        m = len(seq)
        n = int(len(res)/m)
        return res.reshape(n, m).T
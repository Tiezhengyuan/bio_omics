'''
pLDDT per residue local distance testing
'''

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns
import requests
import os



from .pdb import Pdb

class CalculatePlddt(Pdb):
    def __init__(self, pdb_file:str=None, outdir:str=None):
        super().__init__(pdb_file, outdir)

    def complex_ca_bfactor(self) -> dict:
        '''
        scan all models and chains
        '''
        data = {
            'pdb_id': self.structure_id,
            'pdb_file': self.pdb_file,
            'bfactor': [],
        }
        ca_iter = self.c_atom_iter()
        for model, chain, atoms in ca_iter:
            bfactors = []
            for atom in atoms:
                residue_number = int(atom.full_id[3][1])
                # b-factor is not assay data but determined by refinement process
                bfactor = float(atom.get_bfactor())
                bfactors.append((residue_number, bfactor))
            if bfactors:
                df = pd.DataFrame(bfactors, columns=['residue_number', 'bfactor'])
                stats = self.bfactor_statistics(df)
                item = {
                    'model': model,
                    'chain': chain,
                    'bfactor': df,
                    'stats':stats,
                }
                data['bfactor'].append(item)
        return data

    def bfactor_statistics(self, df:pd.DataFrame, res_start:int=None, res_end:int=None) -> dict:
        """
        Calculate B-factor statistics
            low bfactor: more stable, ordered, or rigid part of structure
            high bfactor: greater aomic displacement and disorder.
        """
        if res_start:
            df=df[(df['residue_number']>=res_start)]
        if res_end:
            df=df[(df['residue_number']<=res_end)]
        num_ca = len(df)
        if num_ca > 1:
            scores = df['bfactor']
            stats = {
                'mean': np.mean(scores),
                'median': np.median(scores),
                'std': np.std(scores),
                'min': np.min(scores),
                'max': np.max(scores),
            }
            return stats
        return {}
        
    def complex_ca_plddt(self) -> dict:
        '''
        scan all models and chains
        '''
        data = {
            'structure_id': self.structure_id,
            'pdb_file': self.pdb_file,
            'pLDDT': [],
        }
        ca_iter = self.c_atom_iter()
        for model, chain, atoms in ca_iter:
            plddt = []
            for atom in atoms:
                residue_number = int(atom.full_id[3][1])
                bfactor = float(atom.get_bfactor())
                plddt.append((residue_number, bfactor))
            if plddt:
                plddt = pd.DataFrame(plddt, columns=['residue_number', 'pLDDT'])
                stats = self.plddt_statistics(plddt)
                item = {
                    'model': model,
                    'chain': chain,
                    'plddt': plddt,
                    'stats':stats,
                }
                data['pLDDT'].append(item)
        return data
        
    def monomer_ca_plddt(self) -> pd.DataFrame:
        '''
        only retrieve the first one or that is monomer
        '''
        plddt = []
        ca_iter = self.c_atom_iter()
        for model, chain, atoms in ca_iter:
            for atom in atoms:
                residue_number = int(atom.full_id[3][1])
                bfactor = float(atom.get_bfactor())
                plddt.append((residue_number, bfactor))
            break
        if plddt:
            plddt = pd.DataFrame(plddt, columns=['residue_number', 'pLDDT'])
            return plddt

    def plddt_statistics(self, plddt:pd.DataFrame, res_start:int=None, res_end:int=None) -> dict:
        """
        Calculate comprehensive pLDDT statistics
            Very high (pLDDT > 90)
            High (90 > pLDDT > 70)
            Low (70 > pLDDT > 50)
            Very low (pLDDT < 50)
        """
        if res_start:
            plddt=plddt[(plddt['residue_number']>=res_start)]
        if res_end:
            plddt=plddt[(plddt['residue_number']<=res_end)]
        num_ca = len(plddt)
        if num_ca > 1:
            scores = plddt['pLDDT']
            stats = {
                'mean_plddt': np.mean(scores),
                'median_plddt': np.median(scores),
                'std_plddt': np.std(scores),
                'min_plddt': np.min(scores),
                'max_plddt': np.max(scores),
                'confidence_categories': {
                    'very_high': len(scores[scores >= 90]),
                    'high': len(scores[(scores >= 70) & (scores < 90)]),
                    'low': len(scores[(scores >= 50) & (scores < 70)]),
                    'very_low': len(scores[scores < 50])
                },
                'percentage_high_confidence': len(scores[scores >= 70]) / num_ca * 100
            }
            return stats
        return {}

    @staticmethod
    def plot_per_residue(plddt_data:dict):
        # Per-residue pLDDT plot
        fig, ax = plt.subplots(figsize=(8,4), layout='tight')
        for key, plddt in plddt_data.items():
            ax.plot(plddt['residue_number'], plddt['pLDDT'], label=key, linewidth=1)
        ax.axhline(y=90, color='green', linestyle='--', alpha=0.7)
        ax.axhline(y=70, color='orange', linestyle='--', alpha=0.7)
        ax.axhline(y=50, color='red', linestyle='--', alpha=0.7)
        ax.set_xlabel('Residue Number')
        ax.set_ylabel('pLDDT Score')
        ax.set_title('Per-residue pLDDT Scores')
        ax.legend()
        ax.grid(True, alpha=0.3)
        return fig, ax

    @staticmethod
    def pie_confidence(plddt):
        # Confidence categories pie chart
        fig, ax = plt.subplots(figsize=(6,4), layout='tight')
        confidence_counts = [
            len(plddt[plddt['pLDDT'] >= 90]),
            len(plddt[(plddt['pLDDT'] >= 70) & (plddt['pLDDT'] < 90)]),
            len(plddt[(plddt['pLDDT'] >= 50) & (plddt['pLDDT'] < 70)]),
            len(plddt[plddt['pLDDT'] < 50])
        ]
        labels = ['Very High', 'Hight', 'Low', 'Very Low']
        ax.pie(confidence_counts, labels=labels, autopct='%1.1f%%', startangle=90)
        ax.set_title('Confidence Categories Distribution')
        return fig, ax
    
    @staticmethod
    def plot_distribution(plddt_data):
        # Histogram
        fig, ax = plt.subplots(figsize=(8,4), layout='tight')
        for plddt in plddt_data:
            ax.hist(plddt['pLDDT'], bins=50, edgecolor='black', alpha=0.5)
        ax.axvline(x=90, color='green', linestyle='--', alpha=0.7)
        ax.axvline(x=70, color='orange', linestyle='--', alpha=0.7)
        ax.axvline(x=50, color='red', linestyle='--', alpha=0.7)
        ax.set_xlabel('pLDDT Score')
        ax.set_ylabel('Frequency')
        ax.set_title('pLDDT Distribution')
        ax.grid(True, alpha=0.3)
        return fig, ax



    
    @staticmethod
    def confidence(plddt_score):
        pool = [(90, 100, 'Very High'), (70, 90, 'High'),
            (50, 70, 'Low'), (0, 50, 'Very Low'),]
        for min_score, max_score, label in pool:
            if min_score <= plddt_score <= max_score:
                return label
            
    @staticmethod
    def boxplot_confidence(plddt_data:dict, figsize:tuple=None):
        confidence_data = []
        for key, plddt in plddt_data.items():
            plddt['confidence'] = plddt['pLDDT'].map(CalculatePlddt.confidence)
            plddt['group'] = key
            confidence_data.append(plddt)
        confidence_data = pd.concat(confidence_data, axis=0)

        #
        if figsize is None:
            figsize = (8,4) 
        fig, ax = plt.subplots(figsize=figsize, layout='tight')
        sns.boxplot(confidence_data, x='confidence', y='pLDDT', hue='group', notch=True, ax=ax)
        ax.set_title('pLDDT Scores by Confidence Category')
        ax.set_xlabel('Model Confidence')
        ax.set_ylabel('pLDDT Score')
        ax.grid(True, alpha=0.3)
        return fig, ax

    @staticmethod
    def barplot_confidence(plddt_data:dict, figsize:tuple=None):
        confidence_data = []
        labels = ['Very Low', 'Low', 'High', 'Very High']
        for key, _plddt in plddt_data.items():
            _plddt['confidence'] = _plddt['pLDDT'].map(CalculatePlddt.confidence)
            counts = _plddt['confidence'].value_counts()
            counts = pd.DataFrame(counts).reset_index()
            counts['cum_percent'] = counts['count'].cumsum()/sum(counts['count'])
            counts['group'] = key
            confidence_data.append(counts)
        confidence_data = pd.concat(confidence_data)

        if figsize is None:
            figsize = (8,4) 
        fig, ax = plt.subplots(figsize=figsize, layout='tight')        
        for confidence in labels[::-1]:
            sub = confidence_data[confidence_data['confidence']==confidence]
            sns.barplot(sub, x='group', y='cum_percent', label=confidence, alpha=.4, ax=ax)
        ax.set_title('pLDDT Scores by Confidence Category')
        ax.set_xlabel('Model Confidence')
        ax.set_ylabel('pLDDT Score')
        ax.grid(True, alpha=0.3)
        return fig, ax


    @staticmethod
    def confidence_percentage(plddt, threshold:int):
        names = ('High', 'Low')
        confidence = {i:0 for i in names}
        for score in plddt['pLDDT']:
            if score >= threshold:
                confidence["High"] += 1
            else:
                confidence["Low"] += 1
        n = len(plddt)
        for k in names:
            confidence[k] *= 100 / n
        return confidence

    @staticmethod
    def dot_confidence(plddt_data, figsize:tuple=None, threshold:int=None):
        if threshold is None:
            threshold = 70
        confidence_data =[]
        for key, plddt in plddt_data.items():
            confidence = CalculatePlddt.confidence_percentage(plddt, threshold)
            counts = pd.DataFrame([confidence])
            counts['structure'] = key
            confidence_data.append(counts)
        confidence_data = pd.concat(confidence_data, axis=0)
        # 
        if figsize is None:
            figsize = (6,4) 
        fig, ax = plt.subplots(figsize=figsize, layout='tight')
        sns.scatterplot(confidence_data, x='Low', y='High', alpha=.4, ax=ax)
        ax.set_title('Model Confidence measured by pLDDT')
        ax.set_xlabel(f'Low Confidence (<{threshold}), %')
        ax.set_ylabel(f'High Confidence (>={threshold}), %')
        for i, row in confidence_data.iterrows():
            ax.text(row['Low'], row['High'], row['structure'])
        return fig, ax 
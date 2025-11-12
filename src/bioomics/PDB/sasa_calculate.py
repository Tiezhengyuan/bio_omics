'''
SASA solvent accessible surface area
“A” (Atom), “R” (Residue), “C” (Chain), “M” (Model), or “S” (Structure).
Note: 
    - example: SasaCalculate(pdb_file).cal_freesasa()
    - althernative approach, example:
        SasaCalculate('./pdb/1A14_H_L_N.pdb').run_freesasa_cmd()
'''
from Bio.PDB.SASA import ShrakeRupley
import freesasa
import os
import numpy as np
import pandas as pd
import subprocess

from ..constants.amino_acids import EXTRA_LONGER_NAMES
from .pdb import Pdb
from .analyze_complex import AnalyzeComplex

class SasaCalculate:

    def __init__(self, pdb_file, prob_radious:float=None):
        self.pdb_file = pdb_file
        self.prob_radious = prob_radious if prob_radious else 1.4
    
    def cal_freesasa(self, method:str=None, hetatm:bool=None, hydrogen:bool=None):
        # default method is LeeRichards
        method = freesasa.ShrakeRupley if method == 'S' else freesasa.LeeRichards
        # heteroatoms: water, ions, ligands, modified residues
        hetatm = False if not hetatm else True
        # include hydrogens or not
        hydrogen = False if not hydrogen else True

        try:
            structure = freesasa.Structure(
                self.pdb_file,
                options = {
                    'hetatm' : hetatm,
                    'hydrogen' : hydrogen,
                    'join-models' : False,
                    'skip-unknown' : False,
                    'halt-at-unknown' : False,
                }
            )
            calc_args = freesasa.Parameters({
                'algorithm': method,
                'probe-radius': self.prob_radious,
                'n-slices': 100,
            })
            result = freesasa.calc(structure, calc_args)
            return self.collect_freesasa_result(result)
        except Exception as e:
            print(f'ERROR: freeSASA on {self.pdb_file} failed. error={e}')
        return None

    def collect_freesasa_result(self, result):
        '''
        result is determined by freesasa.calc()
        '''
        try:
            sasa = []
            res_areas = result.residueAreas()
            for chain_id, chain_sasa in res_areas.items():
                for res_area in chain_sasa.values():
                    _sasa = {
                        'chain_id': chain_id,
                        'res_no': res_area.residueNumber,
                        'res': res_area.residueType,
                        'total': res_area.total,
                        'relative_total': res_area.relativeTotal,
                        'apolar': res_area.apolar,
                        'relative_apolar': res_area.relativeApolar,
                        'polar': res_area.polar,
                        'relative_polar': res_area.relativePolar,
                        'side_chain': res_area.sideChain,
                        'relative_side_chain': res_area.relativeSideChain,
                        'main_chain': res_area.mainChain,
                        'relative_main_chain': res_area.relativeMainChain,
                    }
                    _sasa['aa'] = EXTRA_LONGER_NAMES.get(_sasa['res'], '')
                    sasa.append(_sasa)
            df = pd.DataFrame(sasa)
            df = df.set_index(['chain_id', 'res_no'])
            return df
        except Exception as e:
            print(f'ERROR: Collection of SASA based on {self.pdb_file} failed. error={e}')
        return None

    @staticmethod
    def cal_delta_sasa(unbounded_df:pd.DataFrame, bounded_df:pd.DataFrame):
        '''
        calculate delta SASA
        '''
        try:
            matched = bounded_df.index.intersection(unbounded_df.index)
            delta = bounded_df.loc[matched, ['res', 'aa']]
            delta['value'] = unbounded_df.loc[matched, 'total'] - bounded_df.loc[matched, 'total']
            delta = delta.reset_index()

            # check if tight connection between two chains
            #  at least tightly connected 1 residue
            sig = delta[delta['value'] > 0]
            if len(sig) >= 1:
                return delta
        except Exception as e:
            pass
        return None

    @staticmethod
    def parse_seq(df, cutoff:float=None, max_span:int=None):
        '''
        detect continuous segments given cutoff and <max_span
        '''
        if not isinstance(df, pd.DataFrame):
            return None
        if max_span is None:
            max_span = 10
        # values less than 5 in default
        sig = df[df['value'] >= cutoff]

        # group significant point
        group = AnalyzeComplex.group_significant_points(sig, max_span)

        # convert to seq
        motif = []
        for g in group:
            start, end = g[0][0], g[-1][0]
            sub = df.loc[start:end,:]
            rec = {
                'seq': ''.join(sub['aa']),
                'sig_res': len(g),
                'start': g[0][0],
                'end': g[-1][0],
                'cutoff': cutoff,
                'max_span': max_span,
            }
            if 'pair_aa' in list(sub):
                # seq should be continuous seq, but pair_aa may not be continuous
                rec['pair_aa'] = ''.join(sub['pair_aa'])
            motif.append(rec)
        return motif

########################################
    def cal_total_freesasa(self):
        '''
        calculate SASA using freesasa python module
        '''
        structure = freesasa.Structure(self.pdb_file)
        result = freesasa.calc(structure)
        res = freesasa.classifyResults(result, structure)
        res['total'] = result.totalArea()
        return res

    def cal_total_sasa(self) -> dict:
        """
        Calculate total SASA using biopython
        """
        p = Pdb(self.pdb_file)
        sr = ShrakeRupley(probe_radius=self.prob_radious, n_points=100)
        sr.compute(p.structure, level="S")
        return pj.structure.sasa

    def cal_sasa_dict(self) -> dict:
        """
        Calculate SASA for all residues
        """
        p = Pdb(self.pdb_file)
        sr = ShrakeRupley(probe_radius=self.prob_radious, n_points=100)
        sr.compute(p.structure, level="R")

        sasa_dict = {}
        for model in p.structure:
            for chain in model:
                for residue in chain:
                    res_id = (chain.id, residue.id[1])
                    sasa_dict[res_id] = 0
                    for atom in residue:
                        sasa = result.atomArea(atom.serial_number)
                        sasa_dict[res_id] += sasa
        return sasa_dict

    def cal_sasa_df(self, model_id, chain_id):
        '''
        return data frame
        '''
        p = Pdb(self.pdb_file)
        p.load_structure()
        sr = ShrakeRupley(probe_radius=self.prob_radious, n_points=100)
        sr.compute(p.structure, level="C")
    
        for model in p.structure:
            for chain in model:
                if model.id == model_id and chain.id == chain_id:
                    sasa = []
                    for res in chain:
                        sasa.append((res.id, res.sasa))
                    df = pd.DataFrame(sasa, columns=['res_id', 'sasa'])
                    return df
    
    def run_freesasa_cmd(self, outdir, method=None):
        '''
        Calculate SASA using command line
        freesasa should be installed
        '''
        method = '-S' if method == 'S&R' else '-L'
        pdb_name = os.path.splitext(os.path.basename(self.pdb_file))[0]
        outfile = os.path.join(outdir, pdb_name + '.txt')
        # run command
        cmd = [
            'freesasa',
            self.pdb_file,
            '--output=', outfile,
            '--format=seq',
            '-p', 'self.prob_radious',
            method,
        ]
        subprocess.call(cmd)



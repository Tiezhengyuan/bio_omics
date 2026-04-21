'''
SASA solvent accessible surface area
“A” (Atom), “R” (Residue), “C” (Chain), “M” (Model), or “S” (Structure).
formular:
    delta SASA = total of unbounded SASA - total of bounded SASA
'''
from ..common import ProcessPickle
from .sasa_calculate import SasaCalculate


class SasaProcess:

    def __init__(self, prob_radious:float=None, delta_threshold:float=None,\
        max_span:int=None, verbose:bool=None):
        self.prob_radious = prob_radious
        self.delta_threshold = float(delta_threshold) if delta_threshold else 1
        self.max_span = int(max_span) if max_span else 30
        self.verbose = verbose
    
    def process(self, pdb_files_pool:list, outfile:str=None) -> list:
        data = []
        for bounded_pdb_file, unbounded_pdb_files in pdb_files_pool:
            # result is dict
            # keys is bounded in dict, unbouned in dict, delta_total_sasa
            result = self.pair_freesasa(bounded_pdb_file, unbounded_pdb_files)
            if result is not None:
                if self.verbose:
                    delta = result['delta_total_sasa']
                    print(f"dSASA = {delta}, bounded PDB: {bounded_pdb_file}")
                # parse seq/motif
                for unbounded in result['unbounded']:
                    unbounded['motif'] = SasaCalculate.parse_seq(
                        unbounded['delta'],
                        cutoff=self.delta_threshold,
                        max_span=self.max_span
                    )
                data.append(result)
        if data:
            if outfile:
                ProcessPickle(outfile).save_pickle(data)
                if self.verbose:
                    print(f"Export SASA calculation to {outfile}")
        return data

    def pair_freesasa(self, bounded_pdb_file:str, unbounded_pdb_files:list):
        '''
        calculate detaSASA = SASA1+ SASA2 - SASA12
        '''
        # calcualte SASA of complex
        df = SasaCalculate(bounded_pdb_file, self.prob_radious).cal_freesasa()
        if df is None:
            print(f"ERROR: Bounded-SASA of {bounded_pdb_file} failed.")
            return None
        
        # calcualte SASA separately
        unbounded = []
        delta_total_sasa = 0
        for f1 in unbounded_pdb_files:
            df1 = SasaCalculate(f1, self.prob_radious).cal_freesasa()
            if df1 is None:
                if self.verbose:
                    print(f"Warning: unbounded-SASA of {f1} failed.")
            else:
                # calculate delta SASA
                delta1 = SasaCalculate.cal_delta_sasa(df1, df)
                if delta1 is not None:
                    total1 = sum(delta1['value'])
                    delta_total_sasa += total1
                    unbounded.append({
                        'pdb_file': f1,
                        'chain_ids': list(delta1['chain_id'].unique()),
                        'delta_total_sasa': total1,
                        'sasa': df1,
                        'delta': delta1,
                    })

        # export data when only binding are detected
        if unbounded:        
            result = {
                'delta_total_sasa': delta_total_sasa,
                'bounded': {
                    'pdb_file': bounded_pdb_file,
                    'chain_ids': list(df.index.get_level_values(0).unique()),
                    'sasa': df,
                },
                'unbounded': unbounded,
            }
            return result
        return None    


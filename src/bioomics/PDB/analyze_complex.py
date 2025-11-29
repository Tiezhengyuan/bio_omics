'''
'''
import pandas as pd
import re

class AnalyzeComplex:

    @staticmethod
    def cal_span(a, b):
        '''
        input pickle: freesasa.p, or dis.p
        '''
        if isinstance(a, str):
            a = int(re.findall(r'\d+', a)[0])
        if isinstance(b, str):
            b = int(re.findall(r'\d+', b)[0])
        return b - a

    @staticmethod
    def group_significant_points(sig_df:pd.DataFrame, max_span:int):
        group, curr = [], []
        for ix, res_no in sig_df['res_no'].items():
            if curr:
                span = AnalyzeComplex.cal_span(curr[-1][1], res_no)
                if span <= max_span:
                    curr.append((ix, res_no))
                else:
                    group.append(curr)
                    curr = [(ix, res_no),]
            else:
                curr.append((ix, res_no))
        else:
            if curr:
                group.append(curr)
        return group

    @staticmethod
    def get_spans(sig_df:pd.DataFrame):
        # contact residues
        sig = sig_df['res_no']
        #spans
        span = [AnalyzeComplex.cal_span(a, b) for a, b in zip(sig[:-1], sig[1:])]
        # filter
        span = [i for i in span if i > 1]
        return span

    @staticmethod
    def get_fragments(sig_df:pd.DataFrame):
        sig = list(sig_df['res_no'])

        frag, frag_len = [], 1
        for ix in range(0, len(sig)-1):
            span = AnalyzeComplex.cal_span(sig[ix], sig[ix+1])
            if span > 1:
                frag.append(frag_len)
                frag_len = 1
            else:
                frag_len += 1
        else:
            frag.append(frag_len)
            frag_len = 1
        return [i for i in frag if i>1]






    @staticmethod    
    def print_meta(data):
        '''
        # input pickle: freesasa.p, or dis.p
        # print pairwise chain type and chain id
        '''
        for _data in data:
            print('key1', list(_data))
            k1, k2 = _data['type'].split('-')
            c1, c2 = _data['chains'][:2]
            print(f'{k1} chain {c1[0]} ~ {k2} chain {c2[0]}')
            print('---'*10)
    
    @staticmethod
    def build_dataset(data, min_fragment:int=1):
        '''
        # input pickle: sasa_seq.p, or dis_seq.p
        # retrieve dataset from machine learning
        '''
        ds_data = []
        for antibody_name, antibody_rec, receptor_name, receptor_rec in data:
            
            antibody_seq_list = [i['seq'] for i in antibody_rec if \
                len(i['seq']) >= min_fragment]
            receptor_seq_list = [i['seq'] for i in receptor_rec if \
                len(i['seq']) >= min_fragment]
            # both input and output are existing
            if antibody_seq_list and receptor_seq_list:
                # take antigen as input
                rec = {
                    'input': ' | '.join(receptor_seq_list),
                    'input_label': receptor_name,
                    'output': ' | '.join(antibody_seq_list),
                    'output_label': antibody_name,
                }
                ds_data.append(rec)
        return ds_data

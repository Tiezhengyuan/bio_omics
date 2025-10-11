'''
antibody sequence database, https://absd.pasteur.cloud/
'''
from copy import deepcopy
import os
import re
import numpy as np
import pandas as pd
from Bio import SeqIO

from .process_data import ProcessData

class Absd(ProcessData):

    def __init__(self, data_dir:str, verbose:bool=False):
        super().__init__(data_dir)
        self.verbose = verbose
    # TODO
    def __call__(self):
        for specie, path in self.scan_files('fasta'):
            data = {}
            for record in self.scan_data(path):
                _id = self.parse_id(record)
                if _id not in data:
                    data[_id] = []
                rec = {
                }
                data[_id].append(rec)
            outfile = os.path.join(self.data_dir, f"{sepcie}.json")
            self.save_json(data, outfile)
            if self.verbose:
                print(outfile)

    def get_records(self, specie_name, record_id):
        res = []
        for specie, path in self.scan_files('fasta'):
            if specie == 'specie_name':
                for record in self.scan_data(path):
                    if record_id in str(record.id):
                        res.append(record)
        return res

    def scan_files(self, extension:str):
        file_iter = self.recursive_files(self.data_dir)
        for path in file_iter:
            if path.endswith(extension):
                file_name = os.path.basename(path)
                specie = '_'.join(re.split(r'\s|_|\.', file_name)[:2])
                yield specie, path
    
    def scan_data(self, path):
        with open(path, 'r') as f:
            for record in SeqIO.parse(f, 'fasta'):
                yield record
    
    def parse_id(self, record):
        res = str(record.id)
        return res

    def df_specie(self, specie_name):
        for specie, label_text in self.scan_files('_V_labels.csv'):
            if specie == specie_name:
                df = pd.read_csv(label_text)
                df['seq_len'] = df['seq'].map(len)
                chain_map = {'H': 'Heavy', 'L': 'Lambda', 'K': 'Kappa',}
                df['chain_type'] = df['chain_type'].map(lambda x: chain_map[x])
                df = df.sort_values('chain_type')
                return df
    
    def df_combined(self):
        # combine data
        data, info = [], []
        for specie, label_text in self.scan_files('_V_labels.csv'):
            df = self.df_specie(specie)
            # sequence must contain regions of CDR1, CDR2, CDR3
            cdf = df[df['cdr']==123].copy()
            # update df
            cdf['specie'] = specie
            cdf['label'] = 'human' if specie == 'Homo_sapiens' else 'other'
            data.append(cdf)

            # update info
            _info = cdf['chain_type'].value_counts().to_dict()
            _info['raw_num'] = df.shape[0]
            _info['cdr_num'] = cdf.shape[0]
            _info['specie'] = specie
            info.append(_info)
        data = pd.concat(data)
        data = data.sample(frac=1, random_state=1)
        data = data.sort_values(['specie', 'chain_type'])
        info = pd.DataFrame(info)
        return data, info

    def export_datasets(self):
        outdir = os.path.join(self.data_dir, 'labels')
        self.init_dir(outdir)
        cols = [
            'query_seq', 'query_seq_cdr', 'query_seq_cdr_trim', 'query_seq_cdr_mask',
            'seq', 'seq_cdr', 'seq_cdr_trim', 'seq_cdr_mask',
            'gap_seq', 'gap_seq_cdr', 'gap_seq_cdr_trim', 'gap_seq_cdr_mask',
        ]

        # collect data
        data, info = self.df_combined()
        print(info)
        for col in cols:
            self.export_datasets_chain(data, col, outdir)
            self.export_datasets_all(data, col, outdir)
            self.export_expand_datasets_chain(data, col, outdir)

    def export_datasets_chain(self, data, col, outdir):
        # split data by heavy and light chains
        chain_types = {
            'heavy': ['Heavy'],
            'light': ['Kappa', 'Lambda'],
            'kappa': ['Kappa',],
            'lambda': ['Kappa',],
        }
        for key, chain in chain_types.items():
            chain_data = data[data['chain_type'].isin(chain)]

            # export raw to csv
            res = chain_data[[col, 'label']]
            outfile = os.path.join(outdir, f'raw_{key}_{col}.csv')
            self.to_label_csv(res, outfile)

            # export unqiue only to csv
            g = chain_data.groupby(col)
            res = []
            for text, sub in g:
                labels = list(set(sub['label']))
                # remove duplicate seq detected more than one species
                # keep the first one if duplicates belong to one specie
                if len(labels) == 1:
                    res.append([text, labels[0]])
            res = pd.DataFrame(res)
            outfile = os.path.join(outdir, f'unique_{key}_{col}.csv')
            self.to_label_csv(res, outfile)

    def export_datasets_all(self, data, col, outdir):
        # export raw to csv
        res = data[[col, 'label']]
        outfile = os.path.join(outdir, f'raw_all_{col}.csv')
        self.to_label_csv(res, outfile)

        g = data.groupby(col)
        res = []
        for text, sub in g:
            labels = list(set(sub['label']))
            # remove duplicate seq detected more than one species
            # keep the first one if duplicates belong to one specie
            if len(labels) == 1:
                res.append([text, labels[0]])
        res = pd.DataFrame(res)
        outfile = os.path.join(outdir, f'unique_all_{col}.csv')
        # self.to_label_csv(res, outfile)

    def export_expand_datasets_chain(self, data, col, outdir):
        # split data by heavy and light chains
        chain_types = {
            'heavy': ['Heavy'],
            'light': ['Kappa', 'Lambda'],
            'kappa': ['Kappa',],
            'lambda': ['Kappa',],
        }
        for key, chain in chain_types.items():
            df = data[data['chain_type'].isin(chain)]
            df1 = df[df['label']=='human']
            for frac in (2, 10, 50):
                df2 = df[df['label']=='other']
                df2 = df2.sample(frac=frac, replace=True, random_state=1)
                df3 = pd.concat([df1, df2])
                # export raw to csv
                res = df3[[col, 'label']]
                if self.verbose:
                    print('balance datasets:', col, res['label'].value_counts())
                outfile = os.path.join(outdir, f'frac{frac}_{key}_{col}.csv')
                self.to_label_csv(res, outfile)


'''
'''
import pandas as pd
import os
import json
import re

from ..common.process_pickle import ProcessPickle
from ..common.dir import Dir
from .distance_analyze import DistanceAnalyze
from .sasa_analyze import SasaAnalyze
from ..protein.analyze_seq import AnalyzeSeq


class CollectContacts:
    def __init__(self, dir:str):
        self._dir = dir

    def get_raw_pdb(self):
        for path in Dir(self._dir).recrusive_files():
            if path.endswith('ent.gz'):
                pdb_id = os.path.basename(path)
                pdb_id = re.sub(r'^pdb|.ent.gz$', '', pdb_id).upper()
                yield pdb_id, path
    
    def filter_raw_pdb(self, pdb_ids:list):
        for pdb_id, path in self.get_raw_pdb():
            if pdb_id in pdb_ids:
                yield pdb_id, path

    def get_results(self):
        dir_iter = Dir(self._dir).recursive_dirs()
        for path in dir_iter:
            pdb_id = os.path.basename(path)
            yield pdb_id, path

    def filter_dirs(self, pdb_ids:list):
        for pdb_id, path in self.get_results():
            if pdb_id in pdb_ids:
                yield pdb_id, path

    def filter_files(self, pdb_ids:list, file_name:str):
        dir_iter = self.filter_dirs(pdb_ids)
        for pdb_id, _dir in dir_iter:
            path = os.path.join(_dir, file_name)
            if os.path.isfile(path):
                yield pdb_id, path

    def meta(self, pdb_ids):
        path_iter = self.filter_files(pdb_ids, 'meta.json')
        for pdb_id, path in path_iter:
            with open(path, 'r') as f:
                data = json.load(f)
                meta = []
                for rec in data['chains']:
                    df = pd.DataFrame(rec['chains'])
                    df = df.assign(
                        model = rec['model_id'],
                        pdb_id = pdb_id,
                    )
                    meta.append(df)
                meta = pd.concat(meta, axis=0)
                yield meta

    def pair_type(self, pdb_ids):
        path_iter = self.filter_files(pdb_ids, 'meta.json')
        for pdb_id, path in path_iter:
            with open(path, 'r') as f:
                data = json.load(f)
                pair = pd.DataFrame({
                    'sasa': data['pair_sasa'],
                    'dist': data['pair_dist'],
                })
                yield pair

    def pair_sasa_one(self, pdb_ids, type_name:str, idx:int, target_seq:str, chain_len:int=None):
        '''
        idx should be 0 or 1
        '''
        aln = AnalyzeSeq(target_seq)

        names, sasa = [], []
        path_iter = self.filter_files(pdb_ids, 'freesasa.p')
        for pdb_id, path in path_iter:
            data = ProcessPickle(path).load_pickle()
            p = SasaAnalyze(data)
            _sasa = p.get_chain_data(type_name, chain_len)
            for chain_id in _sasa:
                names.append((pdb_id, chain_id[idx]))
                sub = _sasa[chain_id][idx]
                query_seq = ''.join(sub['aa'])
                adj_val = aln.parse_values(query_seq, sub['value'])
                # print(adj_val)
                sasa.append(pd.Series(adj_val))
        
        sasa = pd.concat(sasa, axis=1)
        sasa.columns = names
        print(sasa.shape)
        return sasa.T

    def pair_sasa(self, pdb_ids, type_name:str=None, chain_len:int=None):
        names, antibody_sasa, antigen_sasa = [], [], []
        path_iter = self.filter_files(pdb_ids, 'freesasa.p')
        for pdb_id, path in path_iter:
            p = SasaAnalyze(path)
            sasa = p.get_chain_data(type_name, chain_len)
            for chain_id in sasa:
                # print(sasa[chain_id])
                names.append((pdb_id, chain_id[0], chain_id[1]))
                antibody, antigen = sasa[chain_id]
                antibody_sasa.append(antibody['value'])
                antigen_sasa.append(antigen['value'])
        
        antibody_sasa = pd.concat(antibody_sasa, axis=1)
        antibody_sasa.columns = names

        antigen_sasa = pd.concat(antigen_sasa, axis=1)
        antigen_sasa.columns = names

        print(antibody_sasa.shape, antigen_sasa.shape)
        return antibody_sasa, antigen_sasa

    def pair_dist_one(self, pdb_ids, type_name:str, idx:int, target_seq:str, chain_len:int=None):
        aln = AnalyzeSeq(target_seq)

        names, dist = [], []
        path_iter = self.filter_files(pdb_ids, 'dist.p')
        for pdb_id, path in path_iter:
            print(path)
            p = DistanceAnalyze({'pickle_file': path})
            _dist = p.get_chain_data(type_name, chain_len)
            break
        #     for chain_id in _dist:
        #         names.append((pdb_id, chain_id[idx]))
        #         sub = _dist[chain_id][idx]
        #         query_seq = ''.join(sub['aa'])
        #         adj_val = aln.parse_values(query_seq, -sub['value'])
        #         # print(adj_val)
        #         dist.append(pd.Series(adj_val))
        
        # dist = pd.concat(dist, axis=1)
        # dist.columns = names
        # dist = dist.T
        # print(dist.shape)
        # return dist

    def pair_dist(self, pdb_ids, type_name:str=None, chain_len:int=None):
        names, antibody_dist, antigen_dist = [], [], []
        path_iter = self.filter_files(pdb_ids, 'dist.p')
        for pdb_id, path in path_iter:
            p = DistanceAnalyze(path)
            dist = p.get_chain_data(type_name, chain_len)
            for chain_id in dist:
                names.append((pdb_id, chain_id[0], chain_id[1]))
                antibody, antigen = dist[chain_id]
                antibody_dist.append(antibody['value'])
                antigen_dist.append(antigen['value'])
        
        antibody_dist = pd.concat(antibody_dist, axis=1)
        antibody_dist.columns = names
        antibody_dist = -antibody_dist

        antigen_dist = pd.concat(antigen_dist, axis=1)
        antigen_dist.columns = names
        antigen_dist = -antigen_dist
        
        print(antibody_dist.shape, antigen_dist.shape)
        return antibody_dist, antigen_dist

    def receptor_sasa(self, pdb_ids):
        names, antigen = [], []
        path_iter = self.filter_files(pdb_ids, 'freesasa.p')
        for pdb_id, path in path_iter:
            p = SasaAnalyze(path)
            sasa = p.get_receptor_data()
            for chain_id in sasa:
                names.append((pdb_id, chain_id[0], chain_id[1]))
                _antigen = sasa[chain_id]
                antigen.append(_antigen['value'].iloc[:self.chain_len])
        antigen = pd.concat(antigen, axis=1)
        antigen = antigen.fillna(0)
        antigen.columns = names
        print(antigen.shape)
        return antigen

    def receptor_dist(self, pdb_ids):
        names, antigen = [], []
        path_iter = self.filter_files(pdb_ids, 'dist.p')
        for pdb_id, path in path_iter:
            p = DistanceAnalyze(path)
            dist = p.get_receptor_data()
            for chain_id in dist:
                names.append((pdb_id, chain_id[0], chain_id[1]))
                _antigen = dist[chain_id]
                antigen.append(_antigen['value'].iloc[:self.chain_len])
        antigen = pd.concat(antigen, axis=1)
        antigen = antigen.fillna(0)
        antigen.columns = names
        antigen = -antigen
        print(antigen.shape)
        return antigen


    @staticmethod
    def match_chains(df1, df2):
        '''
        intersection by column names
        '''
        match = df1.columns.intersection(df2.columns)
        df1_match = df1[match].T
        df2_match = df2[match].T
        print(df1_match.shape, df2_match.shape)
        return df1_match, df2_match
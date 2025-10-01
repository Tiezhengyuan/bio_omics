import os
from unittest import TestCase
from ddt import ddt, data, unpack
from pprint import pprint

from src.bioomics import ImgtGenedb

@ddt
class TestImgtGenedb(TestCase):

    def setUp(self):
        data_dir = '/home/yuan/data/IMGT'
        self.obj = ImgtGenedb(data_dir)
    
    def test_parse_genelist(self):
        df = self.obj.parse_genelist(True)
        g = df.groupby(['Species', 'IMGT/GENE-DB'])
        assert set([len(b) for a,b in g]) == {1, }

        human = df[df['Species']=='Homo_sapiens']
        g = human.groupby('IMGT/GENE-DB')
        d = {i:sub.to_dict(orient='records') for i,sub in g}

        results = d['IGHA1']
        assert len(results) == 1
        result = results[0]
        assert result['Species'] == 'Homo_sapiens'
        assert result['IMGT/GENE-DB'] == 'IGHA1'
        assert result['IMGT and HGNC gene definition'] == \
            'Homo sapiens immunoglobulin heavy constant alpha 1'
        assert result['IMGT/LIGM-DB reference sequence(s) for allele *01'] == 'BK063801'
        assert result['NCBI Gene'] == '3493'
        assert result['chain'] == 'heavy'
        assert result['domain'] == 'c'

        results = d['IGHD']
        assert len(results) == 1
        result = results[0]
        assert result['IMGT/LIGM-DB reference sequence(s) for allele *01'] == \
            'K02882;K02881;K02875;K02877;K02876;K02879;K02878;K02880'
        assert result['NCBI Gene'] == '3495'
        assert result['chain'] == 'heavy'
        assert result['domain'] == 'c'
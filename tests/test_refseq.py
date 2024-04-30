'''
Test class ConnectNCBI
'''
from tests.helper import *
from src.bioomics import Refseq

@ddt
class TestRefseq(TestCase):

    @skip
    @data(
        [False, "AP_000001.1"],
        [True, "A0A0D3MVN2"],
    )
    @unpack
    def test_parse_uniprotkb(self, by_uniprotkb, expect):
        res = Refseq(DIR_DATA).parse_uniprotkb(by_uniprotkb)
        assert expect in res
    
    @skip
    def test_parse_epitope(self):
        infile = '/home/yuan/bio/bio_omics/tests/data/NCBI/refseq/mRNA_Prot/human.10.protein.gpff.gz'
        Refseq(DIR_DATA).parse_epitope(infile)


    def test_process_epitope(self):
        entity_path = os.path.join(DIR_DATA, 'epitope')
        Refseq(DIR_DATA).process_epitope(entity_path)


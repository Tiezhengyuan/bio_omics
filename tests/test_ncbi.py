'''
Test class ConnectNCBI
'''
from tests.helper import *
from src.bioomics import NCBI

@ddt
class TestNCBI(TestCase):

    @data(
        ['vertebrate_mammalian', {'vertebrate_mammalian': os.path.join(DIR_DATA, \
            'NCBI/assembly_summary/vertebrate_mammalian/assembly_summary.txt')}],
        ['wrong', {}],
    )
    @unpack
    def test_download_assembly_summary(self, input, expect):
        _, res = NCBI(DIR_DATA).download_assembly_summary([input,])
        assert res == expect

    @data(
        ['https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14',
            'Homo sapiens', 'GRCh38.p14', os.path.join(DIR_DATA, 'NCBI', 'genome', \
            'Homo sapiens', 'GRCh38.p14'), 16],
    )
    @unpack
    def test_download_genome(self, ftp_path, specie, version, expect_path, expect_files):
        local_path, local_files = NCBI(DIR_DATA, False).download_genome(ftp_path, specie, version)
        assert local_path == expect_path
        assert len(local_files) == expect_files

    # TODO: confirm output
    @skip
    def test_download_gene_data(self):
        local_files = NCBI(DIR_DATA).download_gene_data()
        print(local_files)
    
    # TODO: confirm output
    @skip
    def test_download_pubmed(self):
        res = NCBI(DIR_DATA).download_pubmed()
        print(res)



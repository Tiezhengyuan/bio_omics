'''
Test class 
'''
from ..helper import *
from src.bioomics import Uniprot, ParseDatRecord


@ddt
class TestParseDatRecord(TestCase):

    def setUp(self):
        self.c = {}
        infile = os.path.join(DATA_DIR, 'UniProt', 'uniprot_sprot.dat')
        parser  = Uniprot(infile).dat_iter()
        for record in parser:
            self.c[record.id] = record

    @data(
        ['Q6GZX4', 'protein']
    )
    @unpack
    def test_molecule_type(self, input, expect):
        record = self.c[input]
        result = ParseDatRecord(record).molecule_type()
        assert result == expect

    @data(
        ['Q6GZX4', 'Frog virus 3 (isolate Goorha) (FV-3)']
    )
    @unpack
    def test_organism(self, input, expect):
        record = self.c[input]
        result = ParseDatRecord(record).organism()
        assert result == expect


    @data(
        ['Q6GZX4', False],
        ['Q96RD9', True],
    )
    @unpack
    def test_is_imgt(self, input, expect):
        record = self.c[input]
        result = ParseDatRecord(record).is_imgt()
        assert result == expect

    @data(
        ['Q6GZX4', False]
    )
    @unpack
    def test_is_drug(self, input, expect):
        record = self.c[input]
        result = ParseDatRecord(record).is_drug()
        assert result == expect

    @data(
        ['Q6GZX4', ['Activator', 'Reference proteome', 'Transcription', 'Transcription regulation']]
    )
    @unpack
    def test_keywords(self, input, expect):
        record = self.c[input]
        result = ParseDatRecord(record).keywords()
        assert result == expect
        
    @data(
        ['Q6GZX4', {'EMBL': ['AY548484'], 'RefSeq': ['YP_031579.1'], 'SwissPalm': ['Q6GZX4'], \
            'KEGG': ['vg:2947773'], 'Proteomes': ['UP000008770'], 'GO': ['GO:0046782'], \
            'InterPro': ['IPR007031'], 'Pfam': ['PF04947']}],
    )
    @unpack
    def test_dbxrefs(self, input, expect):
        record = self.c[input]
        result = ParseDatRecord(record).dbxrefs()
        assert result == expect
        
    @data(
        ['Q6GZX4', [{'type': 'CHAIN', 'start': '0', 'end': '256', 'seq': 'MAFSAEDVLKEYDRRRRMEAL'+\
            'LLSLYYPNDRKLLDYKEWSPPRVQVECPKAPVEWNNPPSEKGLIVGHFSGIKYKGEKAQASEVDVNKMCCWVSKFKDAMRRY'+\
            'QGIQTCKIPGKVLSDLDAKIKAYNLTVEGVEGFVRYSRVTKQHVAAFLKELRHSKQYENVNLIHYILTDKRVDIQHLEKDLV'+\
            'KDFKALVESAHRMRQGHMINVKYILYQLLKKHGHGPDGPDILTVKTGSKGVLYDDSFRKIYTDLGWKFTPL', \
            'note': 'Putative transcription factor 001R'}]],
    )
    @unpack
    def test_features(self, input, expect):
        record = self.c[input]
        result = ParseDatRecord(record).features()
        assert result == expect

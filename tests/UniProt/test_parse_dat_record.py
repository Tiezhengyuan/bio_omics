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
        # NB: molecule_type is a property (no parentheses)
        record = self.c[input]
        result = ParseDatRecord(record).molecule_type
        assert result == expect

    @data(
        ['Q6GZX4', 'Frog virus 3 (isolate Goorha) (FV-3)']
    )
    @unpack
    def test_organism(self, input, expect):
        record = self.c[input]
        result = ParseDatRecord(record).organism
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

    def test_features_multi(self):
        # Q6GZX3: 4 features, note + evidence qualifiers are flattened
        # to scalars (Biopython stores qualifier values as lists)
        record = self.c['Q6GZX3']
        parser = ParseDatRecord(record)
        seq = parser.seq
        expect = [
            {'type': 'CHAIN', 'start': '0', 'end': '320', 'seq': seq[0:320],
                'note': 'Uncharacterized protein 002L'},
            {'type': 'TRANSMEM', 'start': '300', 'end': '318', 'seq': seq[300:318],
                'note': 'Helical', 'evidence': 'ECO:0000255'},
            {'type': 'REGION', 'start': '260', 'end': '294', 'seq': seq[260:294],
                'note': 'Disordered', 'evidence': 'ECO:0000256|SAM:MobiDB-lite'},
            {'type': 'COMPBIAS', 'start': '261', 'end': '294', 'seq': seq[261:294],
                'note': 'Pro residues', 'evidence': 'ECO:0000256|SAM:MobiDB-lite'},
        ]
        assert parser.features() == expect

    @data(
        ['Q6GZX4', '001R_FRG3G'],
    )
    @unpack
    def test_name(self, input, expect):
        record = self.c[input]
        assert ParseDatRecord(record).name == expect

    @data(
        ['Q6GZX4', 'Q6GZX4'],
        ['Q96RD9', 'Q96RD9'],
    )
    @unpack
    def test_uniprot_acc(self, input, expect):
        record = self.c[input]
        assert ParseDatRecord(record).uniprot_acc == expect

    def test_seq(self):
        res = ParseDatRecord(self.c['Q6GZX4']).seq
        assert res.startswith('MAFSAEDVLKEYDRRRRMEALL')
        assert len(res) == 256
        assert res == res.upper()

    def test_creation_date(self):
        res = ParseDatRecord(self.c['Q6GZX4']).creation_date
        assert res == datetime(2011, 6, 28)

    def test_seq_date(self):
        res = ParseDatRecord(self.c['Q6GZX4']).seq_date
        assert res == datetime(2004, 7, 19)

    def test_dates_missing(self):
        # annotations without dates -> None (no KeyError)
        class FakeRecord:
            annotations = {}
        assert ParseDatRecord(FakeRecord()).creation_date is None
        assert ParseDatRecord(FakeRecord()).seq_date is None

    def test_basic(self):
        res = ParseDatRecord(self.c['Q6GZX4']).basic()
        assert set(res.keys()) == {'id', 'name', 'pro_seq', 'molecule_type',
            'organism', 'dbxrefs', 'features', 'keywords'}
        assert res['id'] == 'Q6GZX4'
        assert res['name'] == '001R_FRG3G'
        assert res['molecule_type'] == 'protein'
        assert res['pro_seq'] == ParseDatRecord(self.c['Q6GZX4']).seq
        # features are exported with flattened qualifiers
        assert res['features'][0]['note'] == 'Putative transcription factor 001R'
        assert res['dbxrefs']['EMBL'] == ['AY548484']

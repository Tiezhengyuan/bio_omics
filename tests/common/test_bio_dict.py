'''
Test src/bioomics/common/bio_dict.py

Ground truth:
- Q6GZX4/Q6GZX3 come from tests/data/UniProt/uniprot_sprot.dat.gz
- the GenBank record is a minimal fixture validated against
  Bio.GenBank.parse (accession is a list, features are Record.Feature
  objects with string location and Qualifier(key, value) pairs)
'''
from ..helper import *
import tempfile
from io import StringIO

from Bio import SeqIO

from src.bioomics import BioDict, BioHandler, Uniprot

GBK_TEXT = '''LOCUS       SCU49845     5028 bp    DNA     circular BCT 30-JUL-2019
DEFINITION  Saccharomyces cerevisiae TPM1 gene.
ACCESSION   SCU49845
VERSION     SCU49845.1
KEYWORDS    .
SOURCE      Saccharomyces cerevisiae
  ORGANISM  Saccharomyces cerevisiae
            Eukaryota; Fungi; Ascomycota.
COMMENT     test comment.
FEATURES             Location/Qualifiers
     source          1..5028
                     /organism="Saccharomyces cerevisiae"
     Region          10..25
                     /note="epitope test"
                     /pseudo
ORIGIN
        1 atgaaattta aaaaatgcaa ttaaaaaaaaa
//'''


class TestBioDict(TestCase):

    @classmethod
    def setUpClass(cls):
        # NB: Q6GZX4 only in the gzipped dat file
        records = {r.id: r for r in
                   Uniprot(os.path.join(DATA_DIR, 'UniProt')).dat_iter()}
        cls.q6gzx4 = records['Q6GZX4']
        # GBK fixture on disk (BioHandler expects a path)
        tmp = tempfile.mkdtemp()
        cls.gbk_file = os.path.join(tmp, 'scu49845.gbk')
        with open(cls.gbk_file, 'w') as f:
            f.write(GBK_TEXT)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(os.path.dirname(cls.gbk_file), ignore_errors=True)

    # ----- swiss_* -----
    def test_swiss_source(self):
        source = BioDict.swiss_source(self.q6gzx4)
        self.assertEqual(source['accession'], 'Q6GZX4')
        self.assertEqual(len(source['seq']), 256)
        refs = source['annotations']['references']
        self.assertIsInstance(refs, list)
        self.assertIn('authors', refs[0])
        self.assertIn('pubmed_id', refs[0])
        # dbxrefs converted to a dict (EMBL: x -> {'EMBL': x})
        self.assertIsInstance(source['dbxrefs'], dict)
        self.assertTrue(len(source['dbxrefs']) > 0)

    def test_swiss_source_malformed_dbxref(self):
        # dbxrefs without ':' must not crash dict()
        record = deepcopy(self.q6gzx4)
        record.dbxrefs = record.dbxrefs + ['malformed-no-colon']
        source = BioDict.swiss_source(record)
        self.assertNotIn('', source['dbxrefs'])

    def test_swiss_feature(self):
        ft = self.q6gzx4.features[0]
        feat = BioDict.swiss_feature(self.q6gzx4, ft)
        self.assertEqual(
            feat['seq'],
            str(self.q6gzx4.seq[ft.location.start:ft.location.end]))
        self.assertLess(feat['seq_start'], feat['seq_end'])
        self.assertIn('qualifiers', feat)
        self.assertEqual(feat['id'], ft.id)

    # ----- gbk_* -----
    def test_gbk_source(self):
        # regression: Bio.GenBank.Record has 'definition', not
        # 'record_definition' (AttributeError)
        record = next(BioHandler.parse_gbk(self.gbk_file))
        source = BioDict.gbk_source(record)
        self.assertEqual(source['locus'], 'SCU49845')
        self.assertEqual(source['accessions'], ['SCU49845'])
        self.assertEqual(
            source['record_definition'],
            'Saccharomyces cerevisiae TPM1 gene')
        self.assertEqual(source['organism'], 'Saccharomyces cerevisiae')
        self.assertTrue(source['seq'].startswith('ATG'))
        self.assertIn('db_source', source)
        self.assertIn('keywords', source)

    def test_gbk_feature(self):
        record = next(BioHandler.parse_gbk(self.gbk_file))
        # Region feature: string qualifier + valueless /pseudo qualifier
        ft = record.features[-1]
        feat = BioDict.gbk_feature(record, ft)
        self.assertEqual(feat['feature_name'], 'Region')
        self.assertEqual(feat['location'], '10..25')
        self.assertEqual(feat['qualifiers']['note'], 'epitope test')
        self.assertEqual(feat['qualifiers']['pseudo'], '')
        self.assertEqual(feat['seq_start'], 10)
        self.assertEqual(feat['seq_end'], 25)
        self.assertEqual(feat['seq'], record.sequence[9:25])

    # ----- fasta -----
    def test_fasta(self):
        record = next(SeqIO.parse(
            StringIO('>acc1.2 some description\nMKVLWA\n'), 'fasta'))
        out = BioDict.fasta(record)
        self.assertEqual(out['id'], 'acc1.2')
        self.assertEqual(out['accession'], 'acc1')
        self.assertEqual(out['name'], 'acc1.2')
        self.assertEqual(out['description'], 'acc1.2 some description')
        self.assertEqual(out['seq'], 'MKVLWA')

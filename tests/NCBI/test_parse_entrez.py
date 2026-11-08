'''
Unit tests for src/bioomics/NCBI/parse_entrez.py (ParseEntrez) and the
offline parts of src/bioomics/NCBI/process_entrez.py (ProcessEntrez).

Records are built in memory (no network access). The db_source annotation
mirrors the real GenPept DBSOURCE line of UniProtKB-derived records, e.g.
    DBSOURCE    UniProtKB: locus APOE_HUMAN, accession P02649; class: standard.
                ...
                xrefs (non-sequence databases): CCDS:CCDS12647.1, PDBsum:1B68,
                PDBsum:1BZ4, ..., Ensembl:ENST00000252486.9, ..., GeneID:348, ...
'''
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

from ..helper import *

from src.bioomics.NCBI.parse_entrez import ParseEntrez
from src.bioomics.NCBI.process_entrez import ProcessEntrez

DB_SOURCE = (
    'UniProtKB: locus APOE_HUMAN, accession P02649; class: standard. '
    'extra accessions:B2RC15,C0JYY5,Q9P2S4 created: Jul 21, 1986. '
    'xrefs: M12529.1, AAB59518.1, '
    'xrefs (non-sequence databases): CCDS:CCDS12647.1, PDBsum:1B68, '
    'PDBsum:1BZ4, PDBsum:1LPE, Ensembl:ENST00000252486.9, '
    'Ensembl:ENSP00000252486.3, GeneID:348, Pfam:PF01242, '
)

PRO_SEQ = (
    'MKVLWAALLVTFLAGCQAKVEQAVETEPEPELRQQTEWQSGQRWELALGRFWDYLRWVQT'
    'LSEQVQEELLSSQVTQELRALMDETMKELKAYKSELEEQLTPVAEETRARLSKELQAAQA'
)


def make_record(record_id='P02649.1'):
    record = SeqRecord(
        Seq(PRO_SEQ),
        id=record_id,
        description='apolipoprotein E precursor',
        annotations={
            'organism': 'Homo sapiens',
            'db_source': DB_SOURCE,
        },
    )
    record.features = [
        SeqFeature(
            FeatureLocation(0, 18),
            type='signal peptide',
            qualifiers={'note': ['Peptide signal']},
        ),
        SeqFeature(
            FeatureLocation(18, 60),
            type='region',
            qualifiers={'region_name': ['Apolipoprotein'], 'note': ['pfam01442']},
        ),
    ]
    return record


class TestParseEntrez(TestCase):

    def setUp(self):
        self.parser = ParseEntrez(make_record())

    def test_get_db_source(self):
        res = self.parser.get_db_source()
        self.assertEqual(res['UniProtKB'], 'P02649')
        # PDBsum ids are collected in order of appearance
        self.assertEqual(res['pdbsum'], ['1B68', '1BZ4', '1LPE'])
        self.assertEqual(
            res['ensembl'], ['ENST00000252486.9', 'ENSP00000252486.3'])
        self.assertEqual(res['geneid'], ['348'])
        self.assertEqual(res['pfam'], ['PF01242'])

    def test_get_db_source_missing_annotation(self):
        record = make_record()
        del record.annotations['db_source']
        self.assertIsNone(ParseEntrez(record).get_db_source())

    def test_get_db_source_refseq_only(self):
        record = make_record()
        record.annotations['db_source'] = 'REFSEQ: accession NM_000041.4'
        res = ParseEntrez(record).get_db_source()
        self.assertEqual(res, {})

    def test_get_features(self):
        res = self.parser.get_features()
        self.assertEqual(len(res), 2)
        signal, region = res
        self.assertEqual(signal['type'], 'signal peptide')
        self.assertEqual(signal['start'], 0)
        self.assertEqual(signal['end'], 18)
        self.assertEqual(signal['seq'], PRO_SEQ[:18])
        self.assertEqual(signal['note'], 'Peptide signal')
        self.assertEqual(region['type'], 'region')
        self.assertEqual(region['seq'], PRO_SEQ[18:60])
        self.assertEqual(region['region_name'], 'Apolipoprotein')
        self.assertEqual(region['note'], 'pfam01442')

    def test_get_features_multi_part_location(self):
        record = make_record()
        record.features = [
            SeqFeature(
                FeatureLocation(0, 3) + FeatureLocation(6, 9),
                type='CDS',
                qualifiers={},
            ),
        ]

        res = ParseEntrez(record).get_features()
        self.assertEqual(len(res), 2)
        self.assertEqual(res[0]['seq'], PRO_SEQ[0:3])
        self.assertEqual(res[1]['seq'], PRO_SEQ[6:9])


def make_entrez_data():
    return {
        'P02649.1': {
            'id': 'P02649.1',
            'desc': 'apolipoprotein E precursor',
            'pro_seq': PRO_SEQ,
            'organism': 'Homo sapiens',
            'db_source': {
                'UniProtKB': 'P02649',
                'pdbsum': ['1BZ4', '1B68'],
            },
            'features': [
                {
                    'type': 'region',
                    'start': 18,
                    'end': 60,
                    'seq': PRO_SEQ[18:60],
                    'region_name': 'Apolipoprotein',
                },
                {
                    'type': 'signal peptide',
                    'start': 0,
                    'end': 18,
                    'seq': PRO_SEQ[:18],
                    'note': 'Peptide signal',
                },
            ],
        },
        'NP_000041.1': {
            'id': 'NP_000041.1',
            'desc': 'apolipoprotein E isoform b',
            'pro_seq': PRO_SEQ,
            'organism': 'Homo sapiens',
            'db_source': {},
            'features': [],
        },
    }


class TestProcessEntrez(TestCase):

    def setUp(self):
        self.pe = ProcessEntrez(make_entrez_data())

    def test_init_defaults(self):
        # no data: empty store
        pe = ProcessEntrez()
        self.assertEqual(pe.data, {})
        self.assertFalse(pe.verbose)
        # data passed in is stored and returned by get_data()
        pe = ProcessEntrez(make_entrez_data())
        self.assertEqual(pe.get_data(), make_entrez_data())

    def test_get_seq_dedups_by_sequence(self):
        seq_id, id_seq = self.pe.get_seq()
        # both records share the same sequence
        self.assertEqual(seq_id[PRO_SEQ], ['P02649', 'NP_000041.1'])
        # UniProtKB accession replaces the record id when present
        self.assertEqual(id_seq['P02649'], PRO_SEQ)
        self.assertEqual(id_seq['NP_000041.1'], PRO_SEQ)

    def test_get_pdb_ids_sorted_and_deduped(self):
        res = self.pe.get_pdb_ids()
        self.assertEqual(res, ['1B68', '1BZ4'])

    def test_get_info_by_uniprot(self):
        rec = self.pe.get_info('P02649')
        self.assertEqual(rec['id'], 'P02649.1')

    def test_get_info_by_record_id(self):
        rec = self.pe.get_info('NP_000041.1')
        self.assertEqual(rec['id'], 'NP_000041.1')

    def test_get_info_missing(self):
        self.assertIsNone(self.pe.get_info('Q99999'))

    def test_cluster_features(self):
        ft_df = self.pe.cluster_features('P02649', ['Apolipoprotein'])
        self.assertIn('ft', ft_df)
        self.assertIn('Apolipoprotein', ft_df)
        self.assertEqual(len(ft_df['ft']), 2)
        self.assertEqual(len(ft_df['Apolipoprotein']), 1)
        self.assertEqual(
            ft_df['Apolipoprotein'].iloc[0]['region_name'], 'Apolipoprotein')
        # the split frames are cached on the record
        rec = self.pe.get_info('P02649')
        self.assertIn('ft_df', rec)
        self.assertEqual(len(rec['ft_df']['ft']), 2)

    def test_cluster_features_unknown_id(self):
        self.assertIsNone(self.pe.cluster_features('Q99999', ['X']))

    def test_cluster_features_no_region_name(self):
        pe = ProcessEntrez({
            'X1': {
                'id': 'X1',
                'pro_seq': 'AAA',
                'organism': 'Homo sapiens',
                'db_source': {},
                'features': [
                    {'type': 'CDS', 'start': 0, 'end': 3, 'seq': 'AAA'}],
            },
        })
        with self.assertRaises(KeyError):
            pe.cluster_features('X1', ['X'])


if __name__ == '__main__':
    import unittest
    unittest.main()


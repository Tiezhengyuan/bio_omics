'''
Test class 
'''
from tests.helper import *

from src.bioomics import IEDBEpitope

@ddt
class TestIEDBEpitope(TestCase):

    # process() requires a ProteinMeta result dir plus the full IEDB
    # epitope/MHC dumps (multi-GB downloads) - network-heavy, skipped
    @skip
    def test_(self):
        IEDBEpitope(DIR_DATA).process()

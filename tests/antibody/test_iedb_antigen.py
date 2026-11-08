'''
Test class 
'''
from tests.helper import *

from src.bioomics import IEDBAntigen

@ddt
class TestIEDBAntigen(TestCase):

    # process() integrates the full IEDB epitope/MHC/B-cell dumps
    # (multi-GB downloads) - network-heavy, skipped
    @skip
    def test_(self):
        IEDBAntigen(DIR_DATA).process()

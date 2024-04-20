'''
Test class 
'''
from tests.helper import *

from src.bioomics import IEDBEpitope

@ddt
class TestIEDBEpitope(TestCase):

    def test_(self):
        IEDBEpitope(DIR_DATA).process()

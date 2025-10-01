'''
Test class 
'''
from tests.helper import *

from src.bioomics import IEDBAntigen

@ddt
class TestIEDBAntigen(TestCase):

    def test_(self):
        IEDBAntigen(DIR_DATA).process()

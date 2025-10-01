'''
Test class 
'''
from .helper import *
from src.bioomics import Expasy


@ddt
class TestExpasy(TestCase):

    def setUp(self):
        self.c = Expasy(DIR_DATA, False)

    def list_files(self):
        res = self.c.list_files()
        print(res)

    def test_download_data(self):
        res = self.c.download_data()
        print(res)

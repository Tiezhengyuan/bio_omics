from ..helper import *
from src.bioomics import ParsePdbLine

print(DATA_DIR)

@ddt
class TestParsePdbLine(TestCase):

    @data(
        ['ATOM      1  N   SER A 606      81.660  -9.273  16.918  1.00 36.35           N', ],
        ['ATOM      2  CA  SER A 606      81.305  -8.249  17.940  1.00 36.46           C ', ],
    )
    @unpack
    def test_parse_atom(self, line):
        result = ParsePdbLine(line).parse_atom()


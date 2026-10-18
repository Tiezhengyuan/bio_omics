'''
example: pytest -s tests/protein/test_plot_seq_region.py
'''
import matplotlib
matplotlib.use('Agg')

from contextlib import redirect_stdout
import io
import matplotlib.pyplot as plt

from ..helper import *

from src.bioomics import PlotSeqRegion

def make_ft():
    return pd.DataFrame([
        {'code': 'W', 'start': 1, 'end': 98, 'length': 98},
        {'code': 'FR1-IMGT', 'start': 1, 'end': 25, 'length': 98},
        {'code': 'CDR1-IMGT', 'start': 26, 'end': 33, 'length': 98},
        {'code': 'FR2-IMGT', 'start': 34, 'end': 50, 'length': 98},
        {'code': 'CDR2-IMGT', 'start': 51, 'end': 58, 'length': 98},
        {'code': 'FR3-IMGT', 'start': 59, 'end': 96, 'length': 98},
    ])

@ddt
class TestPlotSeqRegion(TestCase):

    def test_chain_len_from_ft(self):
        fig, ax = plt.subplots()
        p = PlotSeqRegion(ax, make_ft())
        # chain length is read from the W row (regression: label-based [0])
        assert p.chain_len == 98
        plt.close(fig)

    def test_chain_len_given(self):
        fig, ax = plt.subplots()
        p = PlotSeqRegion(ax, make_ft(), 120)
        assert p.chain_len == 120
        plt.close(fig)

    def test_add_dia(self):
        fig, ax = plt.subplots()
        p = PlotSeqRegion(ax, make_ft(), 98)
        assert p.add_dia('FR1-IMGT') is None
        assert p.lines is not None
        assert p.boxes is not None
        # one box for the region
        assert len(p.boxes) == 1
        # a head line (start > 0) and a tail line (chain_len > end)
        assert len(p.lines) == 2
        assert p.ylim == 2
        plt.close(fig)

    def test_add_dia_wrong_name(self):
        fig, ax = plt.subplots()
        p = PlotSeqRegion(ax, make_ft(), 98)
        buf = io.StringIO()
        with redirect_stdout(buf):
            assert p.add_dia('NOPE') is None
        assert 'wrong region name' in buf.getvalue()
        plt.close(fig)

    def test_set_ax(self):
        fig, ax = plt.subplots()
        p = PlotSeqRegion(ax, make_ft(), 98)
        p.set_ax()
        assert ax.get_xlim() == (-60.0, 118.0)
        plt.close(fig)
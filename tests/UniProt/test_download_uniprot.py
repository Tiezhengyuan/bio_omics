'''
Test class DownloadUniprot: FTP of UniProtKB/Swiss-Prot
'''
import tempfile
from pathlib import Path

from ..helper import *
from src.bioomics import DownloadUniprot
from src.bioomics.connector.conn_ftp import ConnFTP

ENDPOINT_COMPLETE = '/pub/databases/uniprot/current_release/knowledgebase/complete'


@ddt
class TestDownloadUniprot(TestCase):

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.outdir = Path(self._tmp.name)

    def test_init(self):
        c = DownloadUniprot(self.outdir, False)
        # target dir <outdir>/UniProt is created
        assert c.local_path == self.outdir / 'UniProt'
        assert c.local_path.is_dir()
        assert c.url == 'ftp.uniprot.org'

    def test_init_str_outdir(self):
        # outdir may be a str: it used to raise TypeError on str / str
        c = DownloadUniprot(str(self.outdir), False)
        assert c.local_path == self.outdir / 'UniProt'

    def test_uniprot_sprot(self):
        c = DownloadUniprot(self.outdir, False)
        with mock.patch.object(ConnFTP, 'download_file',
                side_effect=lambda endpoint, file_name, local_path:
                    str(c.local_path / file_name)) as m:
            res = c.uniprot_sprot()
        # 3 files are downloaded: dat, fasta, xml
        assert [os.path.basename(f) for f in res] == \
            ['uniprot_sprot.dat.gz', 'uniprot_sprot.fasta.gz', 'uniprot_sprot.xml.gz']
        # each file is fetched from the knowledgebase complete directory
        assert [call.kwargs['file_name'] for call in m.call_args_list] == \
            ['uniprot_sprot.dat.gz', 'uniprot_sprot.fasta.gz', 'uniprot_sprot.xml.gz']
        for call in m.call_args_list:
            assert call.kwargs['endpoint'] == ENDPOINT_COMPLETE
            assert call.kwargs['local_path'] == c.local_path

    def test_uniprot_trembl(self):
        c = DownloadUniprot(self.outdir, False)
        with mock.patch.object(ConnFTP, 'download_file',
                side_effect=lambda endpoint, file_name, local_path:
                    str(c.local_path / file_name)) as m:
            res = c.uniprot_trembl()
        assert len(res) == 3
        assert [call.kwargs['endpoint'] for call in m.call_args_list] == \
            [ENDPOINT_COMPLETE] * 3

    def test_uniprot_sprot_skip_existing(self):
        # existing local files are not downloaded again when overwrite=False
        c = DownloadUniprot(self.outdir, False)
        names = ('uniprot_sprot.dat.gz', 'uniprot_sprot.fasta.gz', 'uniprot_sprot.xml.gz')
        for name in names:
            (c.local_path / name).write_text('data')
        with mock.patch.object(ConnFTP, 'download_file') as m:
            res = c.uniprot_sprot()
        m.assert_not_called()
        assert res == [c.local_path / name for name in names]

    @data(
        ['uniref50', '/pub/databases/uniprot/current_release/uniref/uniref50'],
        ['uniref100', '/pub/databases/uniprot/current_release/uniref/uniref100'],
    )
    @unpack
    def test_uniref(self, prefix, endpoint):
        c = DownloadUniprot(self.outdir, False)
        with mock.patch.object(ConnFTP, 'download_file',
                return_value=str(c.local_path / f'{prefix}.fasta.gz')) as m:
            res = getattr(c, prefix)()
        m.assert_called_once_with(
            endpoint=endpoint,
            file_name=f'{prefix}.fasta.gz',
            local_path=c.local_path,
        )
        assert res == c.local_path / f'{prefix}.fasta.gz'

    @data(
        ['uniref50'],
        ['uniref100'],
    )
    @unpack
    def test_uniref_skip_existing(self, prefix):
        # the local fasta.gz (correct name) is detected: no download
        c = DownloadUniprot(self.outdir, False)
        local_file = c.local_path / f'{prefix}.fasta.gz'
        local_file.write_text('data')
        with mock.patch.object(ConnFTP, 'download_file') as m:
            res = getattr(c, prefix)()
        m.assert_not_called()
        assert res == local_file



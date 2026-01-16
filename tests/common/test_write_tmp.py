from ..helper import *
from src.bioomics.common.write_tmp import WriteTmp
from src.bioomics import ProcessJson

@ddt
class TestParsePdb(TestCase):

    def setUp(self):
        pass
    
    def test_files(self):
        tmp_dir = '/home/yuan/test/tmp'
        pdb_id = '3GKI'
        data_dir = '/home/yuan/test/pdb'
        dst_dir = os.path.join(data_dir, pdb_id[:2], pdb_id)
        with WriteTmp(tmp_dir, dst_dir) as writer:
            print(writer.temp_path)
            (writer.temp_path / 'sub').mkdir()
            for i in range(10):
                outfile = os.path.join(writer.temp_path, f'{i}.json')
                ProcessJson.save_json({i:i}, outfile)


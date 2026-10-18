'''
Test performance of concurrency in Python

Unit tests for src/bioomics/common/concurrency.py; the heavy benchmark
(download of gene_refseq_uniprotkb_collab.gz from NCBI FTP) is kept at
the bottom, skipped.
'''
from ..helper import *
from src.bioomics import *
from concurrent.futures import ThreadPoolExecutor

# module level worker functions: multiprocessing requires picklable targets
def _square(x):
    return x * x

def _touch_file(args):
    outdir, x = args
    with open(os.path.join(outdir, f'{x}.done'), 'w') as f:
        f.write(str(x * x))


# ----- benchmark runners (require the NCBI FTP download) -----
def single(entity_path, nrows, nworks, chunk_size):
    c = IntegrateNCBIProtein(DIR_DATA, False, nworks)
    c.integrate = IntegrateData(entity_path)
    c.meta = c.integrate.get_meta({})
    c.index_meta = c.integrate.get_index_meta()

    # donwload and parsing
    gz_file = NCBI(DIR_DATA).download_gene_map('gene_refseq_uniprotkb_collab.gz')
    data_iter = ParseID(nrows).gene_map0(gz_file, '#NCBI_protein_accession')
    c.integrate_uniprotkb(data_iter)

def mt(entity_path, nrows, nworks, chunk_size):
    c = IntegrateNCBIProtein(DIR_DATA, False, nworks)
    c.integrate = IntegrateData(entity_path)
    c.meta = c.integrate.get_meta({})
    c.index_meta = c.integrate.get_index_meta()

    gz_file = NCBI(DIR_DATA).download_gene_map('gene_refseq_uniprotkb_collab.gz')
    data_pool_iter = ParseID(nrows, nworks, chunk_size).gene_map(gz_file, '#NCBI_protein_accession')
    for data_pool in data_pool_iter:
        multi_threads(data_pool, c.integrate_uniprotkb)

def ex(entity_path, nrows, nworks, chunk_size):
    c = IntegrateNCBIProtein(DIR_DATA, False, nworks)
    c.integrate = IntegrateData(entity_path)
    c.meta = c.integrate.get_meta({})
    c.index_meta = c.integrate.get_index_meta()

    gz_file = NCBI(DIR_DATA).download_gene_map('gene_refseq_uniprotkb_collab.gz')
    data_iter = ParseID(nrows, nworks, chunk_size).gene_map2(gz_file, '#NCBI_protein_accession')
    with ThreadPoolExecutor(nworks) as executor:
        executor.map(c.integrate_uniprotkb, data_iter)

def mp(entity_path, nrows, nworks, chunk_size):
    c = IntegrateNCBIProtein(DIR_DATA, False, nworks)
    c.integrate = IntegrateData(entity_path)
    c.meta = c.integrate.get_meta({})
    c.index_meta = c.integrate.get_index_meta()

    gz_file = NCBI(DIR_DATA).download_gene_map('gene_refseq_uniprotkb_collab.gz')
    data_pool_iter = ParseID(nrows, nworks, chunk_size).gene_map(gz_file, '#NCBI_protein_accession')
    for data_pool in data_pool_iter:
        mp_pool(data_pool, c.integrate_uniprotkb)

def wrapper(func, nrows, nworks, chunk_size):
    entity_path = os.path.join(DIR_DATA, 'demo')
    if os.path.isdir(entity_path):
        os.system(f"rm -fr {entity_path}")
    print(nrows, nworks, chunk_size)

    start = datetime.now()
    func(entity_path, nrows, nworks, chunk_size)
    end = datetime.now()
    delta = end-start
    return delta.seconds


class TestConcurrencyUnit(TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)

    def test_multi_threads(self):
        results = []
        # list.append is atomic under the GIL
        multi_threads([1, 2, 3, 4], lambda x: results.append(_square(x)))
        self.assertEqual(sorted(results), [1, 4, 9, 16])

    def test_multi_threads_empty(self):
        multi_threads([], lambda x: None)  # no crash

    def test_multi_processes(self):
        # NB: target functions must be module level for fork/spawn
        multi_processes([(self.tmp, i) for i in range(4)], _touch_file)
        done = sorted(int(f.split('.')[0]) for f in os.listdir(self.tmp))
        self.assertEqual(done, [0, 1, 2, 3])

    def test_multi_processes_empty(self):
        multi_processes([], _touch_file)  # no crash

    def test_mp_pool(self):
        self.assertEqual(sorted(mp_pool([1, 2, 3, 4], _square)), [1, 4, 9, 16])

    def test_mp_pool_generator(self):
        # NB: the old version called len(data_iter): TypeError on generators
        self.assertEqual(sorted(mp_pool((i for i in range(4)), _square)),
                         [0, 1, 4, 9])

    def test_mp_pool_empty(self):
        self.assertEqual(mp_pool([], _square), [])


class TestConcurrency(TestCase):

    def test_thread_pool_executor(self):
        # documented alternative in the original benchmark
        with ThreadPoolExecutor(4) as executor:
            results = list(executor.map(_square, [1, 2, 3]))
        self.assertEqual(results, [1, 4, 9])

    @skip('performance benchmark: downloads 1M rows from NCBI FTP')
    def test_disk_io(self):
        df = pd.DataFrame({
            'nrows': [i*10000 for i in [1,10,100,100,100]],
            'nworks': [4,8,8,12,16],
        })
        chunk_size = 10000
        df = df.assign(chunk_size=chunk_size, single=0, multithreads=0,
            threads_executor=0, multiprocess=0)

        for row in range(0, df.shape[0]):
            nrows, nworks = df.loc[row, 'nrows'], df.loc[row, 'nworks']
            duration = wrapper(single, nrows, nworks, chunk_size)
            df.loc[(df.nrows==nrows) & (df.nworks==nworks), 'single'] = duration
            duration = wrapper(mt, nrows, nworks, chunk_size)
            df.loc[(df.nrows==nrows) & (df.nworks==nworks), 'multithreads'] = duration
            duration = wrapper(ex, nrows, nworks, chunk_size)
            df.loc[(df.nrows==nrows) & (df.nworks==nworks), 'threads_executor'] = duration
            duration = wrapper(mp, nrows, nworks, chunk_size)
            df.loc[(df.nrows==nrows) & (df.nworks==nworks), 'multiprocess'] = duration

        print(df)



'''
Test performance of concurrency in Python
'''
from tests.helper import *
from src.bioomics import *
from concurrent.futures import ThreadPoolExecutor

def single(entity_path, nrows, nworks, chunk_size):
    c = IntegrateNCBIProtein(DIR_DATA, False, nworks)
    c.integrate = IntegrateData(entity_path)
    c.meta = c.integrate.get_meta({})
    c.index_meta = c.integrate.get_index_meta()

    # donwload and parsing
    gz_file = NCBI(DIR_DATA).download_gene_map('gene_refseq_uniprotkb_collab.gz')
    data_iter = ParseID(nrows).gene_map0(gz_file, '#NCBI_protein_accession')
    c.integrate_uniprotkb(data_iter)

# multithread
def mt(entity_path, nrows, nworks, chunk_size):
    c = IntegrateNCBIProtein(DIR_DATA, False, nworks)
    c.integrate = IntegrateData(entity_path)
    c.meta = c.integrate.get_meta({})
    c.index_meta = c.integrate.get_index_meta()

    # donwload and parsing
    gz_file = NCBI(DIR_DATA).download_gene_map('gene_refseq_uniprotkb_collab.gz')
    data_pool_iter = ParseID(nrows, nworks, chunk_size).gene_map(gz_file, '#NCBI_protein_accession')
    for data_pool in data_pool_iter:
        multi_threads(data_pool, c.integrate_uniprotkb)
# multithread: executor
def ex(entity_path, nrows, nworks, chunk_size):
    c = IntegrateNCBIProtein(DIR_DATA, False, nworks)
    c.integrate = IntegrateData(entity_path)
    c.meta = c.integrate.get_meta({})
    c.index_meta = c.integrate.get_index_meta()

    # donwload and parsing
    gz_file = NCBI(DIR_DATA).download_gene_map('gene_refseq_uniprotkb_collab.gz')
    data_iter = ParseID(nrows, nworks, chunk_size).gene_map2(gz_file, '#NCBI_protein_accession')
    with ThreadPoolExecutor(nworks) as executor:
        executor.map(c.integrate_uniprotkb, data_iter)

# multiprocess
def mp(entity_path, nrows, nworks, chunk_size):
    c = IntegrateNCBIProtein(DIR_DATA, False, nworks)
    c.integrate = IntegrateData(entity_path)
    c.meta = c.integrate.get_meta({})
    c.index_meta = c.integrate.get_index_meta()

    # donwload and parsing
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

class TestConcurrency(TestCase):

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


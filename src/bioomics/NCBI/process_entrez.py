import io
import pandas as pd
from Bio import Entrez, SeqIO

from .parse_entrez import ParseEntrez


Entrez.email = 'tiezhengyuan@hotmail.com'

class ProcessEntrez:
    def __init__(self, data=None, verbose=False):
        self.data = data
        self.verbose = verbose

    def search_protein_ids(self, keyword):
        """
        Search NCBI's protein database for IDs matching a keyword
        """
        handle = Entrez.esearch(
            db="protein",
            term=keyword,
            retmax=1000
        )
        res = Entrez.read(handle)
        handle.close()
        if self.verbose:
            print(f"Found {res['Count']} protein IDs.")
        return res['IdList']

    def fetch_records(self, protein_ids):
        self.data = {}
        for i, _id in enumerate(protein_ids):
            record = self.fetch_record(_id)
            entrez_parser = ParseEntrez(record)
            rec_dict = {
                'id':record.id,
                'desc':record.description,
                'pro_seq': str(record.seq),
                'organism': record.annotations['organism'],
                'db_source': entrez_parser.get_db_source(),
                'features': entrez_parser.get_features(),
            }
            self.data[record.id] = rec_dict
        return self.data

    def fetch_record(self, query_id):
        """
        Fetch complete protein records for given IDs
        """
        handle = Entrez.efetch(
            db="protein",
            id=query_id,
            rettype="gb",
            retmode="text"
        )
        record_text = handle.read()
        handle.close()
        # 
        handle = io.StringIO(record_text)
        record = SeqIO.read(handle, 'genbank')
        return record

    def get_seq(self):
        seq_id, id_seq = {}, {}
        for _id, rec in self.data.items():
            seq = rec['pro_seq']
            if 'db_source' in rec and 'UniProtKB' in rec['db_source']:
                _id = rec['db_source']['UniProtKB']
            if seq not in seq_id:
                seq_id[seq] = []
            if _id not in seq_id[seq]:
                seq_id[seq].append(_id)
            id_seq[_id] = seq
        return seq_id, id_seq

    def get_pdb_ids(self):
        res = []
        for _id, rec in self.data.items():
            if 'db_source' in rec and 'pdbsum' in rec['db_source']:
                pdb = rec['db_source']['pdbsum']
                res.extend(pdb)
        res = list(set(res))
        res.sort()
        return res

    def get_info(self, query_id:str):
        for _id, rec in self.data.items():
            if 'db_source' in rec and 'UniProtKB' in rec['db_source'] \
                and query_id == rec['db_source']['UniProtKB']:
                return rec
            if query_id == _id:
                return rec
        return None

    def cluster_features(self, pdb_id:str, keywords:list):
        '''
        update self.data['ft_df']
        '''
        record = self.get_info(pdb_id)
        ft = pd.DataFrame(record['features'])
        print(record['id'], ft['region_name'].value_counts())

        # split by features
        record['ft_df'] = {
            'ft': ft,
        }
        for k in keywords:
            sub = ft[ft['region_name'].str.contains(k, na=False)]
            record['ft_df'][k] = sub
        return record['ft_df']

    def get_data(self):
        return self.data
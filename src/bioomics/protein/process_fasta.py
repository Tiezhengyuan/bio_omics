from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class ProcessFasta:
    def __init__(self, fa_file:str, verbose:bool=False):
        self.fa_file = fa_file
        self.verbose = verbose

    def one_seq(self, seq:str, seq_id:str):
        rec = SeqRecord(
            Seq(seq),
            id= str(seq_id),
            description=''
        )
        return self.to_fasta([rec,])

    def row_seq(self, rows):
        records = []
        for row in rows:
            rec = SeqRecord(
                Seq(row['seq']),
                id= str(row['id']),
                description=row.get('description', '')
            )
            records.append(rec)
        return self.to_fasta(records)
    
    def unique_sequences(self, df):
        records = []
        g = df.groupby('seq')
        for seq, sub in g:
            rec = SeqRecord(
                Seq(seq),
                id = sub['id'].str.cat(sep=','),
                description=''
            )
            records.append(rec)
        return self.to_fasta(records)

    def to_fasta(self, records):
        try:
            with open(self.fa_file, 'w') as f:
                SeqIO.write(records, f, 'fasta')
            if self.verbose:
                print(f"Export to {self.fa_file}.")
            return 'succeed'
        except Exception as e:
            if self.verbose:
                print(f"Failure: Can't write to {self.fa_file}, error={e}")
        return 'fail'
   
    def load_fasta(self) -> list:
        try:
            with open(self.fa_file, 'r') as f:
                parser = SeqIO.parse(f, 'fasta')
                records = [rec for rec in parser]
                return records
        except Exception as e:
            if self.verbose:
                print(f"Failure: Can't read FASTA {self.fa_file}, error={e}")
        return []

    def get_first_seq(self) -> str:
        for rec in self.load_fasta():
            return str(rec.seq)
        return ''

    def wrong_fasta(self) -> bool:
        '''
        check if a fasta file is good
        false: sequence is empty string, or not at least one record
        '''
        data = self.load_fasta()
        if not data:
            return True
        for record in data:
            seq = str(record.seq)
            if len(seq) == 0:
                return True
        return False

    def clean_fasta(self):
        '''
        remove empty sequence record
        '''
        remove_empty = False
        records = []
        data = self.load_fasta()
        for record in data:
            seq = str(record.seq)
            if len(seq) > 0:
                records.append(record)
            else:
                remove_empty = True
        self.to_fasta(records)
        return remove_empty

from Bio import SeqIO


class ProcessFasta:
    def __init__(self, fa_file:str, verbose:bool=False):
        self.fa_file = fa_file
        self.verbose = verbose

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

    def load_fasta(self):
        try:
            with open(self.fa_file, 'r') as f:
                parser = SeqIO.parse(f, 'fasta')
                records = [rec for rec in parser]
                return records
        except Exception as e:
            if self.verbose:
                print(f"Failure: Can't write to {self.fa_file}, error={e}")
        return None

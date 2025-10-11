from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re

class ProcessSeq:

    @staticmethod
    def insert_chr(seq:str, idx_chr:dict):
        '''
        example indx = {idx: chr}
        '''
        pool, offset = list(seq), 0
        for idx, chr in idx_chr.items():
            pool.insert(idx + offset, chr)
            offset += 1
        return ''.join(pool)
    
    @staticmethod
    def encompass(seq:str, idx_pool:list, left_sep:str, right_sep:str):
        idx_chr = {}
        for start, end in idx_pool:
            idx_chr[start] = left_sep
            idx_chr[end] = right_sep
        seq = ProcessSeq.insert_chr(seq, idx_chr)
        return seq

    @staticmethod
    def encompass_substr(seq:str, substr_pool:list, left_sep:str, right_sep:str):
        idx_chr = {}
        for substr in substr_pool:
            match = re.search(substr, seq)
            if match:
                start, end = match.start(), match.end()
                idx_chr[start] = left_sep
                idx_chr[end]= right_sep
        seq = ProcessSeq.insert_chr(seq, idx_chr)
        return seq

    @staticmethod
    def to_fasta(id2seq:list, outfile:str=None):
        records = []
        try:
            for item in id2seq:
                record = SeqRecord(
                    Seq(item['seq']),
                    id = item['id'],
                    description = item.get('description', '')
                )
                records.append(record)
            if outfile:
                with open(outfile, 'w') as f:
                    SeqIO.write(records, f, 'fasta')
                return True
        except Exception as e:
            print(e)
        return False

    @staticmethod
    def mask_seq(seq:str, regexp:str=None, labels:list=None):
        if regexp is None:
            regexp = r'\[.*?\]'
        if labels is None:
            labels = [f'<CDR{i}>' for i in range(1, 6)]
        else:
            labels = [str(i) for i in labels]
        offset = 0
        mask_seq = seq
        for match, mask in zip(re.finditer(regexp, seq), labels):
            start, end = match.start(), match.end()
            mask_seq = mask_seq[:start+offset] + mask + mask_seq[end+offset:]
            offset -= (end-start)
            offset += len(mask)
        if mask_seq != seq:
            return mask_seq

    @staticmethod
    def digit_cdr_mask(seq_cdr_mask:str):
        if not seq_cdr_mask:
            return seq_cdr_mask
        seq = re.sub(r'(<CDR)(\d)(>)?', r'\2', seq_cdr_mask)
        seq = re.sub(r'\d$', '', seq)
        return seq

    @staticmethod
    def trim_cdr3(seq_cdr:str):
        if not seq_cdr:
            return seq_cdr
        matches = [m for m in re.finditer(r'\[.*?\]', seq_cdr)]
        if len(matches) > 2:
            match_cdr3 = matches[2]
            _end = match_cdr3.start()
            return seq_cdr[:_end]
        return seq_cdr
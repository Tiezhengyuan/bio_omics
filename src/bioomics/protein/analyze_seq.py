import re
import numpy as np
import pandas as pd
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner

class AnalyzeSeq:
    def __init__(self, target_seq):
        self.target_seq = target_seq
        self.aln = []
        self.aligner = PairwiseAligner()
        self.aligner.mode = 'global'  # or 'local'
        self.aligner.match_score = 1
        self.aligner.mismatch_score = 0
        self.aligner.open_gap_score = 0
        self.aligner.extend_gap_score = 0

    def iter_align(self, seq_iter):
        for item in seq_iter:
            alignments = self.align_pairwise(item['query_seq'])
            score = round(alignments.score)
            # exact match: ratio=1
            n1, n2 = len(self.target_seq), len(item['query_seq'])
            ratio = round(score*2/(n1+n2), 4)
            # only consider the top hit in alignments
            rec = {
                'name': item.get('name'),
                'score': score,
                'ratio': ratio,
                'alignment': alignments[0],
            }
            if 'query_val' in item:
                _val = item['query_val']
                rec['query_val'] = _val
                rec['adj_val'] = self.adjust_value(alignments[0], _val)
            self.aln.append(rec)
        return self.aln
    
    def align_pairwise(self, query_seq:str):
        seq1 = Seq(self.target_seq)
        seq2 = Seq(query_seq)
        alignments = self.aligner.align(seq1, seq2)
        return alignments
   
    def adjust_value(self, alignment, val):
        '''
        adjust values of query according to alignment with target
        '''
        target_len = len(alignment.target)
        _min = np.min(list(val))
        adj_val = np.array([_min,] * target_len)
    
        for _target, _query in zip(alignment.aligned[0], alignment.aligned[1]):
            target_indices = np.arange(_target[0], _target[1])
            query_indices = np.arange(_query[0], _query[1])
            adj_val[target_indices] = val[query_indices]
        return adj_val

    def parse_values(self, query_seq, query_val):
        alignments = self.align_pairwise(query_seq)
        adj_val = self.adjust_value(alignments[0], query_val)
        return adj_val

    def print_alignments(self):
        for alignment in alignments:
            print(alignment)

    def parse_igblastp_text(self, output_file):
        results = []
        current = {}
        summary_header = None
        alignment = []
        
        with open(output_file) as f:
            for line in f:
                line = line.strip()
                if not line:
                    if summary_header:
                        alignment = pd.DataFrame(alignment)
                        alignment = alignment.set_index(alignment[0])
                        alignment = alignment.iloc[:,1:]
                        alignment.columns = summary_header
                        current['alignment'] = alignment
                        summary_header = None
                    continue
                if line.startswith("Query="):
                    if current:  # Save previous result
                        results.append(current)
                    # start new query
                    query = line.split()[1].split('|')
                    keys = ['pdb_id', 'chain_id', 'chain_type']
                    current = {a:b for a,b in zip(keys, query)}
                    alignment = []
                elif line.startswith("Length="):
                    current["query_len"] = int(line.split('=')[-1])

                if line.startswith('Alignment summary'):
                    header = re.findall(r'\((.*)\)', line)[0]
                    summary_header = header.split(', ')
                    alignment = []
                if summary_header and '-IMGT' in line:
                    items = line.split('\t')
                    alignment.append(items[:len(summary_header)+1])
        if current:  # Add last result
            results.append(current)
            
        return results

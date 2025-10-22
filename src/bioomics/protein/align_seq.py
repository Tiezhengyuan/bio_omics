from itertools import combinations
import os
import io
import re
import numpy as np
import pandas as pd
import subprocess

from Bio.Seq import Seq


class AlignSeq:

    def __init__(self, verbose:bool=False):
        self.verbose = verbose

    @staticmethod    
    def pairwise_align(seq_dict):
        '''
        globalxx - Simplest global alignment (no parameters, match=+1, mismatch=0)
        '''
        from Bio import pairwise2

        res = []
        seqs = list(seq_dict)
        for _seq1, _seq2 in combinations(seqs, 2):
            seq1, seq2 = Seq(_seq1), Seq(_seq2)
            alignments = pairwise2.align.globalxx(seq1, seq2)
            score = alignments[0].score
            ratio = score*2/(len(_seq1) + len(_seq2))
            rec = {
                'seq1': _seq1,
                'id1': ','.join(seq_dict[_seq1]),
                'seq2': _seq2,
                'id2': ','.join(seq_dict[_seq2]),
                'score': int(score),
                'ratio': round(ratio, 2),
                'alignments': alignments,
            }
            res.append(rec)
        return res
    
    @staticmethod
    def view_alignment(aln):
        from Bio import pairwise2
        for score in aln:
            alignments = score['alignments']
            for alignment in alignments:
                l1, l2, s = len(seq1), len(seq2), alignment.score
                print(f"seq1={l1}, seq2={l2}, score={s}")
                print(pairwise2.format_alignment(*alignment))
                break

    def muscle(self, fa_file:str, outfile:str, large_msa:bool=None):
        # alignment
        method = '-super5' if large_msa else '-align'
        cmd = ['muscle', method, fa_file, '-output', outfile,]
        cmd = ' '.join(cmd)
        try:
            if self.verbose:
                print('MUSCLE for multiple sequence alignment: ', cmd)
            result = subprocess.run(cmd, shell=True, \
                capture_output=True, text=True, check=True)
            if self.verbose:
                print("Command output:", result.stdout.strip())
            return 'succeed'
        except subprocess.CalledProcessError as e:
            if self.verbose:
                print(f"Command failed with error: {e}")
                print(f"Stderr: {e.stderr}")
        return 'fail'
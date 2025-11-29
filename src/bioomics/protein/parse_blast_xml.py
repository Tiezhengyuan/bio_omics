'''
https://biopython.org/docs/dev/api/Bio.Blast.html
'''
import pandas as pd
from Bio.Blast import NCBIXML

class ParseBlastXml:
    def __init__(self, infile):
        # infile is BLAST XML output, outfmt=5
        self.infile = infile

    def best_hits(self):
        results = []
        for record in self.record_iter():
            query = {
                'query': record.query,
                'query_len': record.query_length,
            }
            alignments = self.alignments(record)
            if alignments:
                df = pd.DataFrame(alignments)
                df = df.sort_values('evalue')
                res = df.iloc[0].to_dict()
                res.update(query)
                if res not in results:
                    results.append(res)
            else:
                results.append(query)
        return results

    def record_iter(self):
        try:
            with open(self.infile, 'r') as f:
                parser = NCBIXML.parse(f)
                for record in parser:
                    yield record
        except Exception as e:
            print(f'ParseBlastXml failed: file={self.infile}, error={e}')
            return []
           
    def alignments(self, record):
        results = []
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                aln = {
                    'hit_id': alignment.hit_id,
                    'hit_desc': alignment.hit_def,
                    'alignment_len': alignment.length,
                    'evalue': hsp.expect,
                    'bit_score': hsp.bits,
                    'identity': hsp.identities/hsp.align_length,
                    'align_len': hsp.align_length,
                    'query_seq': hsp.query,
                    'match_seq': hsp.match,
                    'subject_seq': hsp.sbjct,
                }
                results.append(aln)
        return results
            


'''
https://biopython.org/docs/dev/api/Bio.Blast.html
'''
import pandas as pd
from Bio.Blast import NCBIXML

class ParseBlastXml:
    def __init__(self, infile):
        # infile is BLAST XML output, outfmt=5
        self.infile = infile

    def record_iter(self):
        try:
            with open(self.infile, 'r') as f:
                parser = NCBIXML.parse(f)
                for record in parser:
                    yield record
        except Exception as e:
            print(f'ParseBlastXml failed: file={self.infile}, error={e}')
            return []

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

    def hits(self):
        results = []
        for record in self.record_iter():
            query = {
                'query': record.query,
                'query_len': record.query_length,
            }
            alignments = self.alignments(record)
            for aln in alignments:
                aln.update(query)
                results.append(aln)
        return results
           
    def alignments(self, record):
        results = []
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                aln = {
                    'hit_id': alignment.hit_id.strip(),
                    'hit_desc': alignment.hit_def,
                    'alignment_len': alignment.length,
                    'evalue': hsp.expect,
                    'bit_score': hsp.bits,
                    'identities': hsp.identities,
                    'identity': hsp.identities/hsp.align_length,
                    'align_len': hsp.align_length,
                    'query_seq': hsp.query,
                    'query_start': hsp.query_start,
                    'query_end': hsp.query_end,
                    'match_seq': hsp.match,
                    'sbjct_seq': hsp.sbjct,
                    'sbjct_start': hsp.sbjct_start,
                    'sbjct_end': hsp.sbjct_end,
                }
                results.append(aln)
        return results
    
    @staticmethod
    def parse_imgt_allele(input:str=None) -> dict:
        '''
        allele_name, specie, gene_name, gene_family, chain_type, isotype
        '''
        if not input:
            return {}

        res = {}
        items = input.split('|')
        allele_name = items[0]
        res['allele_name'] = allele_name
        if len(items) > 1:
            res['specie'] = items[1]
        if allele_name:
            if '*' in allele_name:
                gene_name = allele_name.split('*')[0]
                res['gene_name'] = gene_name
                if '-' in gene_name:
                    res['gene_family'] = gene_name.split('-')[0]
            res['chain_type'] = allele_name[2]
            res['isotype'] = allele_name[:2]
        return res


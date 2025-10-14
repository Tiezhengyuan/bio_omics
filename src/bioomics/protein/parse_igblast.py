import os
import re
import numpy as np
import numpy as np
import pandas as pd
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
from Bio import SeqIO

from ..common import Utils, ProcessJson, ParseText
from .process_seq import ProcessSeq

class ParseIgblast:

    def __init__(self, aln_file:str, fa_file:str=None, verbose:bool=None):
        self.aln_file = aln_file
        self.fa_file = fa_file
        self.verbose = verbose if verbose else False

    def iter_pairwise(self):
        iter1 = self.iter_query()
        iter2 = self.iter_fmt3()
        for record, block in zip(iter1, iter2):
            yield record, block

    def iter_query(self):
        with open(self.fa_file, 'r') as f:
            parser = SeqIO.parse(f, 'fasta')
            for record in parser:
                yield record

    def iter_fmt3(self):
        with open(self.aln_file, 'r') as f:
            tag, rec = False, []
            for line in f:
                # note: only remove newline but keep white space is there are
                line = line.rstrip('\n')
                if line.startswith('Query='):
                    if rec:
                        yield rec
                    tag = True
                    rec = []
                if tag:
                    rec.append(line)
            if rec:
                yield rec

    def analyze_aln(self, specie=None):
        data = []
        info = {
            'specie': specie,
            'fasta': Utils.count_iterator(self.iter_query()),
            'aln': Utils.count_iterator(self.iter_fmt3()),
            'pairwise': 0,
        }
        n = 0
        if self.verbose:
            print(specie, ':')
        for record, block in self.iter_pairwise():
            info['pairwise'] += 1
            try:
                query_seq = str(record.seq)
                rec = {
                    'record_id': record.id,
                    'query_id': self.query_id(block),
                    'query_length': self.query_length(block),
                    'significant': self.sig_alignment(block),
                }
                alignment = self.alignment_summary(block, query_seq)
                if alignment:
                    rec['alignment'] = alignment
                    # annote CDR based on query seq
                    rec['query'] = self.integrate_cdr(query_seq, alignment)
                    # annote CDR based on alignment seq
                    alignment_seq = self.alignment_seq(query_seq, alignment)
                    rec['alignment'] = self.integrate_cdr(alignment_seq, alignment)
                    # gapped seq based on aligned seq
                    region, aln_seq, pos = self.alignments(block)
                    rec['gap'] = self.gap_seq(region, aln_seq, pos)
                # update data
                data.append(rec)
                n += 1
                if self.verbose and n % 10_000 == 0:
                    print(n, end=',')
            except Exception as e:
                if self.verbose:
                    print(block[0], e)
        
        # save data to .json
        prefix = self.aln_file.replace('.aln', '')
        ProcessJson.save_json(data, f'{prefix}.json')
        info['json'] = len(data)
        if self.verbose:
            print(f'Analyze igblastp result and export to {prefix}.json')
        return info
    
    def build_label_text(self):
        # save sequence to csv
        prefix = self.aln_file.replace('.aln', '')
        input_data = ProcessJson.load_json(f'{prefix}.json')

        data = []
        for rec in input_data:
            cdr = rec['query'].get('cdr', '')
            if rec['significant'] and cdr:
                top = rec['significant'][0]
                query = rec['query']
                align = rec['alignment']
                gap = rec['gap']
                text = {
                    'cdr': cdr,
                    'isotype': top['isotype'],
                    'chain_type': top['chain_type'], 
                    'gene_family': top['gene_family'], 
                    'gene_name': top['gene_name'], 
                    'allele_name': top['allele_name'],
                    # query
                    'query_seq': query['seq'], 
                    'query_seq_cdr': query.get('seq_cdr'),
                    'query_seq_cdr_trim': query.get('seq_cdr_trim'),
                    'query_seq_cdr_mask': ProcessSeq.digit_cdr_mask(query.get('seq_cdr_mask')), 
                    # v-region
                    'seq': align['seq'], 
                    'seq_cdr': align.get('seq_cdr'), 
                    'seq_cdr_trim': align.get('seq_cdr_trim'), 
                    'seq_cdr_mask': ProcessSeq.digit_cdr_mask(align.get('seq_cdr_mask')),
                    # gap seq of v-region
                    'gap_seq': gap['gap_seq'], 
                    'gap_seq_cdr': gap.get('gap_seq_cdr'), 
                    'gap_seq_cdr_trim': gap.get('gap_seq_cdr_trim'), 
                    'gap_seq_cdr_mask': ProcessSeq.digit_cdr_mask(gap.get('gap_seq_cdr_mask')),
                }
                data.append(text)
        data = pd.DataFrame(data)
        outfile = f'{prefix}_labels.csv'
        data.to_csv(outfile, index=False)
        print(f'output: {outfile}')

    def retrieve_hits(self):
        '''
        database is IMGT V-REGION
        prepare records for db insertion
        '''
        hits = []
        for record, block in self.iter_pairwise():
            sig_hits = self.sig_alignment(block)
            if sig_hits:
                sig = sig_hits[0]
                if float(sig['bits-score']) >= 100 and float(sig['e-value']) <= 1e-10:
                    sig['query_id'] = self.query_id(block)
                    sig['query_seq'] = self.query_seq(record)
                    hits.append(sig)
        return hits

    def retrieve_regions(self):
        '''
        database is IMGT V-REGION
        prepare records for db insertion
        '''
        hit_regions = []
        names = ['from', 'to', 'length', 'matches', 'mismatches', 'gaps', 'percent identity']
        for record, block in self.iter_pairwise():
            sig_hits = self.sig_alignment(block)
            if sig_hits:
                sig = sig_hits[0]
                if float(sig['bits-score']) >= 100 and float(sig['e-value']) <= 1e-10:
                    query_id = self.query_id(block)
                    query_seq = self.query_seq(record)
                    alignment = self.alignment_summary(block, query_seq)
                    for region_name, items in alignment.items():
                        rec = {
                            'query_id': query_id,
                            'region': region_name,
                        }
                        rec.update(items)
                        hit_regions.append(rec)
        return hit_regions

    def query_seq(self, record):
        return str(record.seq).upper()

    def query_id(self, block):
        for line in block:
            if line.startswith('Query= '):
                return line.replace('Query= ', '')

    def query_length(self, block):
        for line in block:
            if line.startswith('Length='):
                res = re.findall(r'\d+', line)
                return int(res[0])

    def sig_alignment(self, block):
        tag, sig = 0, []
        for line in block:
            if line.startswith('Sequences producing'):
                tag = 1
                continue
            if line == '':
                tag += 1
                continue
            if tag >=3:
                break
            if tag == 2:
                items = re.split(r'\s{2,}', line)
                keys = ['allele_name', 'bits-score', 'e-value']
                if ' ' in items[0]:
                    items = items[0].split(' ', 1) + items[1:]
                    keys.insert(1, 'product')
                rec = {k:v for k,v in zip(keys, items)}
                rec['gene_name'] = rec['allele_name'].split('*')[0]
                rec['gene_family'] = rec['gene_name'].split('-')[0]
                rec['isotype'] = rec['gene_name'][:2]
                rec['chain_type'] = rec['gene_name'][2]
                sig.append(rec)
        return sig

    def alignment_summary(self, block, query_seq:str=None):
        columns, alignment = None, {}
        for line in block:
            if line.startswith('Alignment summary '):
                header = re.findall(r'\((.+)\)', line)[0]
                columns = header.split(', ')
            else:
                if columns:
                    if line:
                        items = line.split('\t')
                        key = re.findall(r'.*-IMGT', items[0])[0] if 'IMGT' in items[0] else items[0]
                        val = [0 if i=='N/A' else int(float(i)) for i in items[1:]]
                        alignment[key] = {a:b for a,b in zip(columns, val)}
                        try:
                            if key == 'Total':
                                alignment[key]['seq'] = query_seq if query_seq else ''
                            else:
                                start, end = alignment[key]['from'] -1, alignment[key]['to']
                                alignment[key]['seq'] = query_seq[start:end] if query_seq else ''
                        except Exception as e:
                            pass
                    else:
                        break
        return alignment


    def integrate_cdr(self, seq, alignment):
        '''
        build query sequence with CDR labeled or with CDR masked
        '''
        seq_cdr, offset1 = None, 0
        cdr_mask, offset2 = None, 0
        for key, aln in alignment.items():
            start = aln['from']
            end = aln['to']
            if key.startswith('CDR'):
                if seq_cdr is None:
                    seq_cdr = list(seq)
                if cdr_mask is None:
                    cdr_mask = list(seq)
                # label CDR regions by []
                seq_cdr.insert(start-1+offset1, '[')
                offset1 += 1
                seq_cdr.insert(end+offset1, ']')
                offset1 += 1
                # mask CDR regions
                k = re.findall(r'CDR\d', key)[0]
                mask = f'<{k}>'
                cdr_mask = cdr_mask[:start+offset2] + list(mask) + \
                    cdr_mask[end+offset2:]
                offset2 = offset2 + len(mask)- (end -start)
        seq_cdr = ''.join(seq_cdr) if seq_cdr else ''
        seq_cdr_mask = ''.join(cdr_mask) if cdr_mask else ''
        cdr = re.findall(r'<CDR(\d*)?>', seq_cdr_mask)
        seq_cdr_mask = re.sub(r'(<CDR)(\d)(>)', r'\2', seq_cdr_mask)
        return {
            'seq': seq,
            'seq_cdr': seq_cdr,
            'seq_cdr_trim': ProcessSeq.trim_cdr3(seq_cdr),
            'seq_cdr_mask': re.sub(r'\d$', '', seq_cdr_mask),
            'cdr': ''.join(cdr),
        }

    def alignment_seq(self, seq, alignment):
        pos = []
        for key, aln in alignment.items():
            start = aln['from']
            end = aln['to']
            if start > 0:
                pos.append(start)
            if end > 0:
                pos.append(end)
        if pos:
            start, end = min(pos) - 1, max(pos)
            return seq[start:end]
        return ''

    def alignments(self, block):
        '''
        gap_seq with or without CDR labled
        '''
        aln_block = ParseText.split_block(block, 'Alignments')
        # gapped_seq_cdr
        region, aln_seq, pos = [], [], []
        for sub_block in aln_block:
            rec = ParseText.parse_lines_items(sub_block)
            # CDR are detected
            if 'Query_' in sub_block[1]:
                region.append(rec[0][-2])
                aln_seq.append(rec[1][-2])
                pos += [int(rec[1][-3]), int(rec[1][-1])]
            else:
                aln_seq.append(rec[0][-2])
                pos += [int(rec[0][-3]), int(rec[0][-1])]
        return region, aln_seq, pos

    def gap_seq(self, region, aln_seq, pos):
        '''
        get gap seq
        '''
        if not pos:
            return {}

        res = {
            'start': min(pos) -1,
            'end': max(pos),
        }
        region = ''.join(region)
        # print('=', region, '=')
        gap_seq = ''.join(aln_seq)
        if gap_seq:
            res['gap_seq'] = gap_seq
            res['seq'] = gap_seq.replace('-', '')
            if region and re.findall(r'<C.*>', region):
                idx_chr = {}
                for match in re.finditer(r'<C.*?>', region):
                    start, end = match.start(), match.end()
                    idx_chr[start] = '['
                    idx_chr[end] = ']'
                if region.endswith(' '):
                    start, end = len(region.rstrip()), len(region)
                    idx_chr[start] = '['
                    idx_chr[end] = ']'
                # gap_seq_cdr
                gap_seq_cdr = ProcessSeq.insert_chr(gap_seq, idx_chr)
                res['gap_seq_cdr'] = ''.join(gap_seq_cdr).strip()
                res['seq_cdr'] = res['gap_seq_cdr'].replace('-', '')
                # cdr mask
                cdr_mask = ProcessSeq.mask_seq(res['gap_seq_cdr'], None, range(1, 6))
                res['gap_seq_cdr_mask'] = re.sub(r'\d$', '', cdr_mask)
                # seq_cdr by trimming CDR3
                res['gap_seq_cdr_trim'] = ProcessSeq.trim_cdr3(res['gap_seq_cdr'])
        return res


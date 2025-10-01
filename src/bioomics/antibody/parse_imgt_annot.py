import os
import re
import gzip

class ParseImgtAnnot:

    def __init__(self, infile):
        self.infile = infile
        self.lines = self._load()

    def _load(self):
        if os.path.isfile(self.infile):
            with gzip.open(self.infile, 'rt') as f:
                lines = f.read().split('\n')
                return lines

    def __call__(self):
        # the order of methods does matter
        res = {
            'inn_file': self.infile,
            'inn_number': self.next_line('INN number'),
            'common_name': self.next_line('Common name'),
            'inn_name': self.next_line('INN name'),
            'receptor_type': self.next_line('IMGT receptor type'),
            'receptor_description': self.next_line('IMGT receptor description'),
            'species': self.next_line('Species'),
            'proposed_list': self.next_line('Proposed list'),
            'recommended_list': self.next_line('Recommended list'),
            'cas_number': self.next_line('Cas number'),
            'molecular_formula': self.next_line('Molecular formula'),
            'glycosylation_sites': self.next_line('Glycosylation sites'),
        }
        # chains
        res['chains'] = []
        for chain_id in self.next_line('Chain ID', []):
            chain_block = self.chain_block(chain_id)
            chain = {
                'chain_id': chain_id,
                'chain_description': self.line_2cols('IMGT chain description', chain_block),
            }
            chain.update(self.chain_sequence(chain_block))
            # v-domain
            vdomain_blocks = self.chain_vdomain_block(chain_block)
            chain['v_domain'] = self.chain_v_domain(vdomain_blocks)
            # c-domain
            cdomain_blocks = self.chain_cdomain_block(chain_block)
            chain['c_domain'] = self.chain_c_domain(cdomain_blocks)
            res['chains'].append(chain)
        # disulfide_bridges
        res['disulfide_bridges'] = {
            'intrachain': self.block_lines('Intrachain disulfide bridges'),
            'interchain': self.block_lines('Interchain disulfide bridges'),
        }
        return res

    def next_line(self, key, default=None):
        for i, line in enumerate(self.lines):
            if key in line:
                _line = self.lines[i+1].strip().replace('REMARK 410 ', '')
                values = _line.split(',')
                return [i.strip() for i in values]
        return default

    def line_2cols(self, key, lines=None):
        if lines is None:
            lines = self.lines
        for line in lines:
            if key in line:
                start = line.index(key)
                s = line[start+len(key):].strip()
                if s:
                    return s
    
    def chain_block(self, chain_id):
        if not chain_id:
            return []
        flag, block = False, []
        for line in self.lines:
            if 'Chain ID' in line and chain_id in line:
                flag = True
            if flag and line == 'REMARK 410':
                break
            if flag:
                _line = line.replace('REMARK 410 ', '')
                block.append(_line)
        return block

    def chain_sequence(self, lines=None):
        if lines is None:
            lines = self.lines
        # collect block lines
        tag, block1, block2 = False, [], []
        for line in lines:
            if 'Chain amino acid sequence' in line:
                tag = True
                continue
            if 'V-DOMAIN' in line:
                break
            if tag is True:
                _line = line.strip().replace('REMARK 410 ', '')
                if '[' in _line or ']' in _line:
                    block1.append(_line)
                else:
                    block2.append(_line)

        # parse data
        sequence = {
            'chain_seq': ''.join(block2),
            'regions': [],
        }
        annot = ''.join(block1).replace('][', ']=[').split('=')
        annot = [re.sub(r'\s', '', s) for s in annot]
        annot = [i for i in annot if i]
        for s in annot:
            region = {
                'type': re.findall(r'^\[?(.*?)[\(\]]', s)[0],
            }
            start = re.findall(r'\((\d*)-', s)
            if start:
                region['start'] = int(start[0])
            end = re.findall(r'-(\d*)\)', s)
            if end:
                region['end'] = int(end[0])
            if start and end:
                start, end = region['start']-1, region['end']-1
                region['seq'] = sequence['chain_seq'][start:end]
            abbr = re.findall(r'\[(\w*)\]\]$', s)
            if abbr:
                region['region'] = abbr[0]
            sequence['regions'].append(region)
        return sequence

    def chain_v_domain(self, vdomain_blocks):
        domains = []
        for block in vdomain_blocks:
            domain = {}
            # split lines
            tag = False
            block_header, block_label, block_seq = [], [], []
            for line in block:
                if 'CDR1' in line:
                    tag = True
                if tag:
                    if '[' in line or ']' in line:
                        block_label.append(line)
                    else:
                        block_seq.append(line)
                else:
                    block_header.append(line)
            # parse header
            for _line in block_header:
                a, b = re.split(r'\s{2,}', _line, maxsplit=1)
                b = b.strip()
                if a not in domain:
                    domain[a] = b
                else:
                    domain[a] += ' ' + b
            if 'IMGT gene and allele' in domain:
                ss = domain['IMGT gene and allele'].replace(') ', '), ')
                domain['IMGT gene and allele'] = ss.split(', ')
            # parse sequence
            if block_seq:
                domain_align = ''.join(block_seq)
                domain['seq_align'] = domain_align
                domain['seq'] = domain_align.replace('.', '')
                domain_span = ''.join(block_label)
                names = re.findall(r'CDR[1-5]', domain_span)
                starts = [i.start() for i in re.finditer(r'\[', domain_span)]
                ends = [i.end() for i in re.finditer(r'\]', domain_span)]
                cdr = []
                for a,b,c in zip(names, starts, ends):
                    cdr.append({
                        'name': a,
                        'start': b,
                        'end': c,
                        'align': domain_align[b:c],
                        'seq': domain_align[b:c].replace('.', ''),
                    })
                domain['cdr'] = cdr
                # insert '[' and ']' into gap align seq
                align_cdr, offset1 = list(domain_align), 0
                align_cdr_mask, offset2 = list(domain_align), 0
                for a, b,c in zip(names, starts, ends):
                    # encompass CDR sequence with [,]
                    align_cdr.insert(b+offset1, '[')
                    offset1 += 1
                    align_cdr.insert(c+offset1, ']')
                    offset1 += 1
                    # mask CDR sequence
                    mask = f"<{a}>"
                    align_cdr_mask.insert(b+offset2, mask)
                    offset2 += 1
                    align_cdr_mask = align_cdr_mask[:b+offset2] + \
                        align_cdr_mask[c+offset2:]
                    offset2 -= (c-b)
                domain['seq_align_cdr'] = ''.join(align_cdr)
                domain['seq_cdr'] = domain['seq_align_cdr'].replace('.', '')
                domain['seq_align_cdr_mask'] = ''.join(align_cdr_mask)
                domain['seq_cdr_mask'] = domain['seq_align_cdr_mask'].replace('.', '')
            domains.append(domain)
        return domains

    def chain_vdomain_block(self, chain_block) -> list:
        start, rec, blocks = None, [], []
        for line in chain_block:
            if 'V-DOMAIN' in line:
                if 'IMGT domain description' in line:
                    start = line.index('IMGT domain description')
                    if rec:
                        blocks.append(rec)
                        rec = []
                if start:
                    rec.append(line[start:])
                else:
                    rec.append(line)
            else:
                if rec:
                    blocks.append(rec)
                rec = []
        if rec:
            blocks.append(rec)
            rec = []
        return blocks

    def chain_c_domain(self, cdomain_blocks):
        domains = []
        for block in cdomain_blocks:
            domain = {
                'seq_align': '',
            }
            for line in block:
                if line.startswith('IMGT'):
                    a, b = re.split(r'\s{2,}', line, maxsplit=1)
                    b = b.strip()
                    if a not in domain:
                        domain[a] = b
                    else:
                        domain[a] += ' ' + b
                else:
                    domain['seq_align'] += line
            domain['seq'] = domain['seq_align'].replace('.', '')
            if 'IMGT gene and allele' in domain:
                ss = domain['IMGT gene and allele'].replace(') ', '), ')
                domain['IMGT gene and allele'] = ss.split(', ')
            domains.append(domain)
        return domains

    def chain_cdomain_block(self, chain_block) -> list:
        start, rec, blocks = None, [], []
        for line in chain_block:
            if 'C-DOMAIN' in line:
                if 'IMGT domain description' in line:
                    start = line.index('IMGT domain description')
                    if rec:
                        blocks.append(rec)
                        rec = []
                if start:
                    rec.append(line[start:])
                else:
                    rec.append(line)
        if rec:
            blocks.append(rec)
            rec = []
        return blocks

    def block_lines(self, key):
        if not key:
            return {}
        tag, res = False, {}
        for line in self.lines:
            if key in line:
                tag = True
                continue
            if line == 'REMARK 410':
                tag = False
            if tag is True:
                _line = line.replace('REMARK 410 ', '')
                kv = _line.split(': ', 1)
                if len(kv) >=2:
                    res[kv[0]] = kv[1]
                else:
                    res[kv[0]] = ''
        return res
    



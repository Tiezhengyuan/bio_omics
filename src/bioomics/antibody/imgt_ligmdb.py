import os
import re
from pprint import pprint

from .process_data import ProcessData

class ImgtLigmdb(ProcessData):

    def __init__(self, data_dir:str, json_file:str=None):
        super().__init__(data_dir, json_file)
        self.data_dir = os.path.join(self.data_dir, 'LIGM-DB')
        # keys: specie, gene_name, region_name
        # self.data = {}
    
    def iter_imgt_dat(self, infile=None):
        if infile is None:
            infile = os.path.join(self.data_dir, 'imgt.dat')
        with open(infile, 'r') as f:
            rec, block = [], []
            for line in f:
                line = line.rstrip()
                if line == '//':
                    if block:
                        rec.append(block)
                    yield rec
                    rec, block = [], []
                elif line == 'XX':
                    rec.append(block)
                    block = []
                else:
                    block.append(line)
    
    def parse_record(self, blocks):
        record = {
            'id': self._id(blocks),
            'date': self._date(blocks),
            'description': self._description(blocks),
            'keywords': self._keywords(blocks),
            'features': self._features(blocks),
            'sequence': self._sequence(blocks),
        }
        return record
    
    def _id(self, blocks):
        for block in blocks:
            if block[0].startswith('ID'):
                for line in block:
                    tag, value = line[:2], line[2:].strip()
                    res = value.split('; ')
                    return res

    def _date(self, blocks):
        for block in blocks:
            if block[0].startswith('DT'):
                res = []
                for line in block:
                    tag, value = line[:2], line[2:].strip()
                    res.append(value)
                return res

    def _description(self, blocks):
        for block in blocks:
            if block[0].startswith('DE'):
                res = []
                for line in block:
                    tag, value = line[:2], line[2:].strip()
                    res.append(value)
                return ' '.join(res)

    def _keywords(self, blocks):
        for block in blocks:
            if block[0].startswith('KW'):
                res = []
                for line in block:
                    tag, value = line[:2], line[2:].strip()
                    if value.endswith('.'):
                        value = value[:-1]
                    res.append(value)
                res = ' '.join(res).split('; ')
                # print(res)
                return res

    def _features(self, blocks):
        res = []
        for block in blocks:
            block = [i for i in block if i.startswith('FT')]
            for line in block:
                items = re.split(r'\s{2,}', line)
                if len(items) == 3:
                    res.append(items[1:])
                elif len(items) == 2:
                    val = items[1]
                    if val.startswith('/'):
                        res[-1].append(val[1:])
                    else:
                        res[-1][-1] += ' ' + val
            if res:
                features = []
                for items in res:
                    # print(items)
                    ft_name = items[0]
                    pos = items[1].split('..')
                    ft = {
                        'name': ft_name,
                        'start': pos[0],
                        'end': pos[1] if len(pos)>1 else None
                    }
                    for val in items[2:]:
                        # print(val)
                        if '=' in val:
                            k, v = re.findall(r'(.*)=(.*)', val)[0]
                            v = v.replace('"', '')
                            if k == 'translation':
                                v = v.replace(' ', '')
                            if k not in ft:
                                ft[k] = v
                            else:
                                if isinstance(ft[k], list):
                                    ft[k].append(v)
                                else:
                                    ft[k] = [ft[k], v]
                        else:
                            ft[val] = True
                    # print(ft)
                    features.append(ft)
                return features

    def _sequence(self, blocks):
        res = {}
        for block in blocks:
            if block[0].startswith('SQ'):
                header = block[0][2:].strip()
                res['header'] = header.split('; ')
                seq = ''
                for line in block[1:]:
                    line = line.strip()
                    seq += re.split(r'\s{2,}', line)[0]
                res['seq'] = seq.replace(' ', '').upper()
                return res
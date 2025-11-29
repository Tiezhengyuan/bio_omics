'''
'''
import pandas as pd


class ParseText:
    
    @staticmethod
    def split_block(block:list, header:str):
        tag = 0
        res, rec = [], []
        for i, line in enumerate(block):
            if line == header:
                tag = 1
                continue
            if tag:
                if block[i] == '' and block[i+1] == '':
                    break
                else:
                    if line:
                        rec.append(line)
                    else:
                        if rec:
                            res.append(rec)
                            rec = []
        if rec:
            res.append(rec)
            rec = []
        return res


    @staticmethod
    def parse_lines_items(lines:list):
        '''
        1. text columns by 1-many whitespace
        2. missing data are whitespace
        '''
        subdf = pd.DataFrame([list(i) for i in lines])
        match = subdf.apply(lambda x: x.str.isspace().all(), axis=0)
        width = max([len(i) for i in lines])
        indices = list(match[match==True].index) + [width,]
        rec = []
        for line in lines:
            start, items = 0, []
            for end in indices:
                if end>start:
                    val = line[start:end]
                    items.append(val)
                start = end + 1
            rec.append(items)
        return rec
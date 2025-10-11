import re

class ParseEntrez:
    def __init__(self, record):
        self.record = record

    def get_db_source(self):
        if 'db_source' not in self.record.annotations:
            return None

        db_resource = {}
        for sub in self.record.annotations['db_source'].split('; '):
            if ':' in sub:
                a, b = sub.split(': ', 1)
                if a == 'UniProtKB':
                    _b = re.findall(r'accession (.+)', b)
                    if _b:
                        db_resource[a] = _b[0]
                id_names = ('PDBsum', 'Ensembl', 'GeneID', 'Pfam')
                for _name in id_names:
                    _ids = re.findall(f'{_name}:([\w|\.]+), ', b)
                    if _ids:
                        db_resource[_name.lower()] = _ids

        # print(record.id, db_resource)
        return db_resource

    def get_features(self):
        pro_seq = str(self.record.seq)
        res = []
        for ft in self.record.features:
            for locus in ft.location.parts:
                ft_rec = {
                    'type': ft.type,
                    'start': int(locus.start),
                    'end': int(locus.end),
                    'seq': pro_seq[int(locus.start):int(locus.end)],
                }
                for k in ('region_name', 'note'):
                    if k in ft.qualifiers:
                        ft_rec[k]= '|'.join(ft.qualifiers[k])
                res.append(ft_rec)
        return res

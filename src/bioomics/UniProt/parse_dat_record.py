'''
record is defined in *.dat
'''

class ParseDatRecord:

    def __init__(self, record):
        self.record = record

    def basic(self) -> dict:
        return  {
            'id': self.uniprot_acc,
            'name': self.name,
            'pro_seq': self.seq,
            'molecule_type': self.molecule_type,
            'organism': self.organism,
            'dbxrefs': self.dbxrefs(),
            'features': self.features(),
            'keywords': self.keywords(),
        }
    
    @property
    def uniprot_acc(self):
        return self.record.id

    @property
    def name(self):
        return str(self.record.name)

    @property
    def seq(self):
        return str(self.record.seq).upper()

    @property
    def molecule_type(self):
        return self.record.annotations['molecule_type']

    @property
    def organism(self):
        return self.record.annotations['organism']
    
    def is_imgt(self) -> bool:
        kw = self.keywords()
        for k in kw:
            if 'Immunoglobulin' in k:
                return True
        return False
    
    def is_drug(self) -> bool:
        for k in self.record.dbxrefs:
            if 'drug' in k.lower():
                return True
        return False

    def keywords(self) -> list:
        annot = self.record.annotations
        return annot.get('keywords', [])

    def dbxrefs(self) -> dict:
        res = {}
        for ref in self.record.dbxrefs:
            if ':' in ref:
                k, v = ref.split(':', 1)
                if k not in res:
                    res[k] = [v,]
                else:
                    res[k].append(v)
        return res

    def features(self) -> list:
        res = []
        for ft in self.record.features:
            loc = ft.location
            rec = {
                'type': ft.type,
                'start': str(loc.start),
                'end': str(loc.end),
            }
            # start and end is integer mostly, but sometime is string
            try:
                rec['seq'] = str(loc.extract(self.record.seq))
            except Exception as e:
                pass
            if ft.qualifiers:
                rec.update(ft.qualifiers)
            res.append(rec)
        return res

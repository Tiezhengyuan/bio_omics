

class ParseAbsd:
    def __init__(self, record):
        self.record = record
    

    def absd_id(self):
        return self.record.id.split("|||")[0]

    def db_source(self):
        dxrefs = []
        names= ['header','dxref_id','vgene_family','source',]
        desc = self.record.description.split("|||")
        for _desc in desc[1:]:
            items = _desc.split(";")
            rec = {k:v for k,v in zip(names, items)}
            dxrefs.append(rec)

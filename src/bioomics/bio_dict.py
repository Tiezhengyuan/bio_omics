'''
convert object of biopython to dictionary
'''

class BioDict:

    @staticmethod
    def reference(references:list):
        '''
        Bio.ExPASy.Prodoc.Reference
        '''
        return [{
            'authors': ref.authors,
            'comment': ref.comment,
            'consrtm': ref.consrtm,
            'journal': ref.journal,
            'location': ref.location, 
            'medline_id': ref.medline_id,
            'pubmed_id': ref.pubmed_id,
            'title': ref.title,
            } for ref in references]

    @staticmethod
    def swiss_source(record):
        '''
        '''
        annot = {}
        for k, v in record.annotations.items():
            if k == 'references':
                annot[k] = BioDict.reference(v)
            else:
                annot[k] = v
        refs = dict([tuple(i.split(':', 1)) for i in record.dbxrefs])
        return {
            'accession': record.id,
            'seq': str(record.seq),
            'annotations': annot,
            'dbxrefs': refs,
        }

    @staticmethod
    def feature(record, ft):
        return {
            'id': ft.id,
            'qualifiers': ft.qualifiers,
            'seq_start': ft.location.start,
            'seq_end': ft.location.end,
            'seq': str(record.seq[ft.location.start:ft.location.end]),
        }
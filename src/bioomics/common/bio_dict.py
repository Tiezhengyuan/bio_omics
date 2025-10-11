'''
convert object of biopython to dictionary
'''
import re

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
    def swiss_feature(record, ft):
        return {
            'id': ft.id,
            'qualifiers': ft.qualifiers,
            'seq_start': ft.location.start,
            'seq_end': ft.location.end,
            'seq': str(record.seq[ft.location.start:ft.location.end]),
        }
    
    @staticmethod
    def gbk_source(record):
        '''
        *.GBK
        '''
        return {
            'accessions': record.accession,
            'db_source': record.db_source,
            'record_definition': record.record_definition,
            'keywords': record.keywords,
            'locus': record.locus,
            'organism': record.organism,
            'seq': record.sequence,
        }

    @staticmethod
    def gbk_feature(record, ft):
        '''
        *.GBK
        '''
        qualifiers = dict([(re.sub(r'^/|=$', '', q.key), \
            re.sub(r'"', '', q.value)) for q in ft.qualifiers])
        # print(record.locus, ft, ft.location)
        locations = re.findall('\d+', ft.location)
        start = int(locations[0]) if len(locations) > 0 else 0
        end = int(locations[1]) if len(locations) > 1 else 0
        return {
            'feature_name': ft.key,
            'location': ft.location,
            'qualifiers': qualifiers,
            'seq_start': start,
            'seq_end': end,
            'seq': record.sequence[start-1:end],
        }

    @staticmethod
    def fasta(record):
        output =  {
            'description': record.description,
            'id': record.id,
            'name': record.name,
            'seq': str(record.seq),
        }
        for attr in ['annotations', 'features', 'dbxrefs']:
            if hasattr(record, attr):
                output[attr] = getattr(record, attr)
        output['accession'] = record.id.split('.', 1)[0]
        return output
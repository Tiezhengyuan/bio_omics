'''
records reday for table insertion
'''
import os
from ..common import ProcessPickle
from ..PDB import AnalyzeDssp, Pdb
from ..protein import ParseBlastXml, ParseIgblast
from ..UniProt import ParseDatRecord
from .query_complex import QueryComplex

class BuildRecords:

    def __init__(self, verbose):
        self.verbose = verbose

    def pdb(self, pdb_file):
        p = Pdb(pdb_file, None, self.verbose)
        p.load_structure()
        pdb_id = p.structure_id

        pdb_url = f"https://www.rcsb.org/structure/{pdb_id}"
        release_date = p.structure.header['release_date']
        name = p.structure.header['name']
        rec = (pdb_id, pdb_url, release_date, name)
        return [rec,]

    def pdb_compound(self, args):
        pdb_id = args['pdb_id']
        pdb_file = args['pdb_file']
        # records
        p = Pdb(pdb_file, None, self.verbose)
        p.load_structure()

        records = []
        for compound_no, cp in p.compound.items():
            compound_no = int(compound_no)
            chains = cp.get('chain').replace(' ', '')
            _text = []
            for _key, _value in cp.items():
                if _key not in ('chain',) and _value:
                    _text.append(f"{_key}: {_value}")
            compound_text = ' ; '.join(_text)
            rec = (pdb_id, compound_no, chains, compound_text)
            records.append(rec)
        return records

    def pdb_method(self, args):
        pdb_id = args['pdb_id']
        pdb_file = args['pdb_file']
        p = Pdb(pdb_file, None, self.verbose)
        p.load_structure()
        method = p.structure_method
        resolution = p.resolution
        if method or resolution:
            rec = (pdb_id, method, resolution)
            return [rec,]
        return []
        
    def pdb_keyword(self, args):
        pdb_id = args['pdb_id']
        pdb_file = args['pdb_file']
        p = Pdb(pdb_file, None, self.verbose)
        p.load_structure()
        records = []
        for keyword in p.keywords:
            rec = (pdb_id, keyword)
            if rec not in records:
                records.append(rec)
        return records
            
    def pdb_chain(self, pdb_file):
        p = Pdb(pdb_file, None, self.verbose)
        p.load_structure()
        records = []
        for chain_id, pdb_id, model_no, chain_no, chain in p.parse_chains():
            # the column pro_seq is not null
            chain_seq = p._chain_seq(chain)
            rec = (chain_id, pdb_id, model_no, chain_no, chain_seq) if chain_seq \
                else (chain_id, pdb_id, model_no, chain_no, None)
            records.append(rec)
        if self.verbose:
            print(f"Warning: no chains {pdb_id} are detectin in {pdb_file}")
        return records
    
    def pdb_chainid(self, args):
        '''
        chain_id is unique and could be null
        '''
        pdb_id = args['pdb_id']
        pdb_file = args['pdb_file']
        p = Pdb(pdb_file, None, self.verbose)
        p.load_structure()
        records = []
        for chain_id, pdb_id, model_no, chain_no, chain in p.parse_chains():
            rec = (pdb_id, model_no, chain_no, chain_id)
            records.append(rec)
        if records:
            return records
        if self.verbose:
            print(f"Warning: no chains are detectin in {pdb_id}:{pdb_file}")
        return [(pdb_id, None, None, None),]

    def pdb_chain_seq(self, args):
        pdb_id = args['pdb_id']
        pdb_file = args['pdb_file']
        p = Pdb(pdb_file, None, self.verbose)
        p.load_structure()
        records = []
        for chain_id, _pdb_id, model_no, chain_no, chain in p.parse_chains():
            # the column pro_seq is not null
            chain_seq = p._chain_seq(chain)
            pro_len = len(chain_seq)
            rec = (chain_id, pro_len, chain_seq) if chain_seq \
                else (chain_id, None, None)
            records.append(rec)
        if not records and self.verbose:
            print(f"Warning: no chains {pdb_id} are detectin in {pdb_file}")
        return records

                 
    def pdb_bfactor(self, args):
        records = []
        pdb_id = args['pdb_id']
        pfile = args['pickle_file']
        try:
            data = ProcessPickle(pfile).load_pickle()
            for item in data.get('bfactor', []):
                model_no = item['model'].id
                chain_no = item['chain'].id
                avg = item['stats'].get('mean')
                rec = (pdb_id, model_no, chain_no, avg)
                if rec not in records:
                    records.append(rec)
        except Exception as e:
            if os.path.isfile(pfile):
                os.remove(pfile)
            print(f"Warning: pfile={pfile}, error={e}")
        if records == []:
            records = [(pdb_id, None, None, None,), ]
        return records

    def pdb_secondary(self, args):
        pdb_id = args['pdb_id']
        # note: pdb_file contains only one chain
        pdb_file = args['pdb_file']
        records = []
        try:
            ad = AnalyzeDssp(pdb_file)
            for model, dssp_data in ad.dssp_iter():
                model_no = model.id
                for chain_no, dssp in dssp_data.items():
                    dssp_seq = ''.join([i['dssp_code'] for i in dssp])
                    rec = (pdb_id, model_no, chain_no, dssp_seq)
                    if rec not in records:
                        records.append(rec)
        except Exception as e:
            print(f"ERROR: Fail in pdb_secondary(). error={e}")
        if records == []:
            records = [(pdb_id, args['model_no'], args['chain_no'], None),]
        return records

    def align_uniprot(self, args):
        pdb_id = args['pdb_id']
        aln_file = args['aln_file']
        records = []
        best_hits = ParseBlastXml(aln_file).best_hits()
        for hit in best_hits:
            chain_id = hit['query']
            rec0 = (chain_id, pdb_id)
            if 'hit_desc' in hit:
                uniprot_acc = hit['hit_desc'].split('|')[1]
                rec = rec0 + (
                    uniprot_acc,
                    float(hit['evalue']),
                    int(hit['bit_score']),
                    float(hit['identity'])
                )
                if rec not in records:
                    records.append(rec)
            else:
                # insert null if not hits detected
                rec0 += (None,) * 4
                if rec0 not in records:
                    records.append(rec0)
        return records

    def uniprot(self, record):
        p = ParseDatRecord(record)
        rec = (p.uniprot_acc, p.name, p.organism, len(p.seq), p.seq,)
        return [rec,]
    
    def uniprot_keyword(self, record):
        uniprot_acc = record.id
        records = []
        keywords = ParseDatRecord(record).keywords()
        for keyword in keywords:
            rec = (uniprot_acc, keyword)
            if rec not in records:
                records.append(rec)
        return records

    def uniprot_dbxrefs(self, record):
        uniprot_acc = record.id
        records = []
        dbxrefs = ParseDatRecord(record).dbxrefs()
        for k, values in dbxrefs.items():
            for v in values:
                rec = (uniprot_acc, k, v)
                if rec not in records:
                    records.append(rec)
        return records
    
    def uniprot_feature(self, record):
        records = []
        uniprot_acc = record.id
        features = ParseDatRecord(record).features()
        for feature in features:
            rec = (
                uniprot_acc,
                feature.get('type'),
                feature.get('start'),
                feature.get('end'),
                feature.get('seq'),
                feature.get('note'),
                feature.get('evidence'),
            )
            records.append(rec)
        if not records:
            rec = (uniprot_acc,) + (None,) * 6
            records.append(rec)
        return records
    
    def uniprot_alphafold_plddt(self, pLDDT_data):
        uniprot_acc = pLDDT_data['uniprot_acc']
        pfile = pLDDT_data['pfile']
        records = []
        for item in pLDDT_data['pLDDT']:
            model_no = item['model'].id
            chain_no = item['chain'].id
            avg_plddt = item['stats'].get('mean_plddt')
            rec = (uniprot_acc, model_no, chain_no, avg_plddt, pfile)
            records.append(rec)
        return records

    def align_imgt(self, args:dict):
        pdb_id = args['pdb_id']
        faa_file = args['faa_file']
        aln_file = args['aln_file']
        pi = ParseIgblast(aln_file, faa_file, self.verbose)
        records = []
        for hit in pi.retrieve_hits():
            chain_id = hit['query_id']
            rec0 = (chain_id, pdb_id)
            if 'allele_name' in hit:
                rec = rec0 + (
                    hit['allele_name'].split('|')[0],
                    hit['gene_name'],
                    hit['gene_family'],
                    hit['chain_type'].upper(),
                    hit['isotype'].upper(), 
                    float(hit['e-value']),
                    float(hit['bits-score']),
                )
                if rec not in records:
                    records.append(rec)
            else:
                # insert null if not hits detected
                rec0 += (None,) * 7
                if rec0 not in records:
                    records.append(rec0)
        return records

    def align_vregion(self, args:dict):
        pdb_id = args['pdb_id']
        faa_file = args['faa_file']
        aln_file = args['aln_file']
    
        pi = ParseIgblast(aln_file, faa_file, self.verbose)
        records = []
        for region in pi.retrieve_regions():
            chain_id = region['query_id']
            pdb_id = chain_id.split('_')[0]
            if 'percent identity' in region:
                rec = (
                    chain_id,
                    pdb_id,
                    region['region'],
                    region['percent identity'],
                    region['from'],
                    region['to'],
                    region['seq'],
                )
                if rec not in records:
                    records.append(rec)
            else:
                # region name is composite primary key
                rec = (chain_id, pdb_id, 'UN', None, None, None, None)
                if rec not in records:
                    records.append(rec)
        return records

    def imgt_genelist(self, record:dict):
        names = ['Species', 'IMGT/GENE-DB', 'IMGT gene functionality', \
            'IMGT and HGNC gene definition', 'Number of alleles', \
            'Chromosome', 'Chromosomal localization', \
            'IMGT/LIGM-DB reference sequence(s) for allele *01', \
            'NCBI Gene', 'chain', 'domain']
        rec = [record.get(i, '') for i in names]
        return [rec,]

    def imgt_inn(self, inn_data):
        names = ('inn_number', 'cas_number', 'common_name', 'inn_name',
            'proposed_list', 'recommended_list', 'receptor_type',
            'receptor_description', 'species',
            'molecular_formula', 'glycosylation_sites')
        rec = [inn_data['inn_file'], ]
        for name in names:
            val = inn_data.get(name)
            if isinstance(val, list):
                rec.append(','.join(val))
            else:
                rec.append(None)
        return [rec,]
    
    def imgt_inn_chain(self, inn_data):
        inn_number = ','.join(inn_data['inn_number'])
        records = []
        for chain in inn_data.get('chains', []):
            seq = chain.get('chain_seq')
            pro_len = len(seq) if seq else None
            rec = (
                inn_number,
                chain['chain_id'],
                chain['chain_description'],
                pro_len,
                seq,
            )
            records.append(rec)
        return records
    
    def imgt_inn_region(self, inn_data):
        records = []
        for chain in inn_data.get('chains', []):
            for region in chain.get('regions', []):
                rec = (
                    chain['chain_id'],
                    region.get('region'),
                    region.get('type'),
                    region.get('start'),
                    region.get('end'),
                    region.get('seq'),
                )
                records.append(rec)
        return records
    
    def imgt_inn_cdomain(self, inn_data):
        records = []
        for chain in inn_data.get('chains', []):
            for domain in chain.get('c_domain', []):
                rec = (
                    chain['chain_id'],
                    domain.get("IMGT domain description"),
                    ','.join(domain.get("IMGT gene and allele", [])),
                    domain.get('seq'),
                    domain.get('seq_align'),
                )
                records.append(rec)
        return records

    def imgt_inn_vdomain(self, inn_data):
        records = []
        for chain in inn_data.get('chains', []):
            for domain in chain.get('v_domain', []):
                rec = (
                    chain['chain_id'],
                    domain.get("IMGT domain description"),
                    domain.get("CDR-IMGT length"),
                    domain.get("FR-IMGT length"),
                    ','.join(domain.get("IMGT gene and allele", [])),
                    domain.get('seq'),
                    domain.get('seq_cdr'),
                    domain.get('seq_align'),
                    domain.get('seq_align_cdr'),
                    domain.get('seq_cdr_mask'),
                    domain.get('seq_align_cdr_mask'),
                )
                records.append(rec)
        return records
    
    def imgt_inn_cdr(self, inn_data):
        records = []
        for chain in inn_data.get('chains', []):
            for domain in chain.get('v_domain', []):
                for cdr in domain.get('cdr', []):
                    start, end = cdr.get("start"), cdr.get("end")
                    len = None
                    try:
                        len = end - start + 1
                    except Exception as e:
                        pass
                    rec = (
                        chain['chain_id'],
                        cdr.get("name"),
                        cdr.get("seq"),
                        cdr.get("align"),
                        start,
                        end,
                        len,
                    )
                    records.append(rec)
        return records
    
    def imgt_inn_blastp(self, aln_file):
        # parse pro_ids
        qc = QueryComplex(self.verbose)
        query_pro_ids = qc.parse_pro_id('imgt_inn_chain')
        ref_pro_ids = qc.parse_pro_id('pdb_chain_seq')
        
        records = []
        best_hits = ParseBlastXml(aln_file).best_hits()
        for hit in best_hits:
            query_id = int(hit['query'])
            for inn_chain_id in query_pro_ids[query_id]:
                rec0 = (query_id, inn_chain_id)
                if 'hit_desc' in hit:
                    ref_id = int(hit['hit_desc'])
                    for pdb_chain_id in ref_pro_ids[ref_id]:
                        rec = rec0 + (
                            ref_id,
                            pdb_chain_id,
                            float(hit['evalue']),
                            int(hit['bit_score']),
                            float(hit['identity'])
                        )
                        records.append(rec)
                else:
                    # insert null if not hits detected
                    rec0 += (None,) * 5
                    records.append(rec0)
        return records

    def ab_absd_source(self, seq_record):
        # 'header','dxref_id','vgene_family','source'
        records = []
        desc = seq_record.description.split("|||")
        for _desc in desc[1:]:
            items = [seq_record.id, desc[0]] + _desc.split(";")
            records.append(items)
        return records

    def absd_align_imgt(self, args:dict):
        seq_id = args['seq_id']
        faa_file = args['faa_file']
        aln_file = args['aln_file']
        pi = ParseIgblast(aln_file, faa_file, self.verbose)
        records = []
        for hit in pi.retrieve_hits():
            if 'allele_name' in hit:
                rec = (
                    seq_id,
                    hit['allele_name'].split('|')[0],
                    hit['gene_name'],
                    hit['gene_family'],
                    hit['chain_type'].upper(),
                    hit['isotype'].upper(), 
                    float(hit['e-value']),
                    float(hit['bits-score']),
                )
                if rec not in records:
                    records.append(rec)
        if not records:
            # insert null if not hits detected
            rec0 = (seq_id,) + (None,) * 7
            records.append(rec0)
        return records

    def absd_align_vregion(self, args:dict):
        seq_id = args['seq_id']
        faa_file = args['faa_file']
        aln_file = args['aln_file']
    
        pi = ParseIgblast(aln_file, faa_file, self.verbose)
        records = []
        for region in pi.retrieve_regions():
            if 'percent identity' in region:
                rec = (
                    seq_id,
                    region['region'],
                    region['percent identity'],
                    region['from'],
                    region['to'],
                    region['seq'],
                )
                if rec not in records:
                    records.append(rec)
        if not records:
            # region name is composite primary key
            rec = (seq_id, 'UN', None, None, None, None)
            records.append(rec)
        return records

        
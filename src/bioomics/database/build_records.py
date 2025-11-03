'''
records reday for table insertion
'''
import os
from ..common import ProcessPickle
from ..PDB import AnalyzeDssp, Pdb
from ..protein import ParseBlastXml, ParseIgblast

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
            if chain_seq:
                rec = (chain_id, pdb_id, model_no, chain_no, chain_seq)
                records.append(rec)
        if self.verbose:
            print(f"Warning: no chains {pdb_id} are detectin in {pdb_file}")
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
           
    def pdb_bfactor(self, args):
        records = []
        pdb_id = args['pdb_id']
        pfile = args['pfile']
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
        return records

    def pdb_secondary(self, args):
        pdb_id = args['pdb_id']
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
            records = [(pdb_id, None, None, None),]
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
            rec0 = (chain_id, pdb_id,)
            if 'percent identity' in region:
                rec = rec0 + (
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
                rec0 += ('UN', None, None, None, None)
                if rec0 not in records:
                    records.append(rec0)
        return records

    def imgt_genelist(self, record:dict):
        names = ['Species', 'IMGT/GENE-DB', 'IMGT gene functionality', \
            'IMGT and HGNC gene definition', 'Number of alleles', \
            'Chromosome', 'Chromosomal localization', \
            'IMGT/LIGM-DB reference sequence(s) for allele *01', \
            'NCBI Gene', 'chain', 'domain']
        rec = [record.get(i, '') for i in names]
        return [rec,]

    def imgt_gene(self, record:dict):
        names = ['germline_id', 'specie', 'isotype', 'chain_type', \
            'gene_name', 'allele_name', 'region_name']
        rec = [record.get(i, '') for i in names]
        return [rec,]

    def imgt_geneseq(self, record:dict):
        names = ['germline_id', 'seq_type', 'seq']
        rec = [record.get(i, '') for i in names]
        return [rec,]
'''
records reday for table insertion
'''
import os
from ..common import ProcessPickle
from ..PDB import AnalyzeDssp, Pdb

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
        rec = (pdb_id, pdb_url, pdb_file, release_date, name)
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
            
    def pdb_chain(self, args):
        pdb_id = args['pdb_id']
        pdb_file = args['pdb_file']
        p = Pdb(pdb_file, None, self.verbose)
        p.load_structure()
        records = []
        chain_ver = {}
        for model_chains in p.get_chains():
            model_no = int(model_chains['model_no'])
            for chain in model_chains['chains']:
                # chain.id is identical in some PDB.
                chain_no = chain['chain_no']
                chain_id_prefix = f"{pdb_id}_{model_no}_{chain_no}"
                key = chain_id_prefix.upper()
                if key not in chain_ver:
                    chain_ver[key] = 0
                else:
                    chain_ver[key] += 1
                chain_id = f"{chain_id_prefix}_{chain_ver[key]}"
                # the column pro_seq is not null
                chain_seq = chain.get('seq')
                if chain_seq:
                    rec = (chain_id, pdb_id, model_no, chain_no, chain_seq)
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
        ad = AnalyzeDssp(pdb_file)

        records = []
        for model, dssp_data in ad.dssp_iter():
            model_no = model.id
            for chain_no, dssp in dssp_data.items():
                dssp_seq = ''.join([i['dssp_code'] for i in dssp])
                rec = (pdb_id, model_no, chain_no, dssp_seq)
                if rec not in records:
                    records.append(rec)
        return records
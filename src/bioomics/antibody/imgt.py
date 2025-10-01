from Bio import SeqIO
from copy import deepcopy
import os
import re
import requests
import pandas as pd
import subprocess
import sys

from .parse_xml import ParseXml
from .process_data import ProcessData
from .parse_imgt_annot import ParseImgtAnnot
from .process_seq import ProcessSeq

class Imgt(ProcessData):
    url = 'https://www.imgt.org/'

    def __init__(self, data_dir:str):
        super().__init__(data_dir)
    
    '''
    3D-DB: INN
    '''
    def retrieve_inn(self):
        s, f = 0, 0
        data = []
        indir = os.path.join(self.data_dir, '3Dstructure-DB', 'IMGT3DFlatFiles')
        for gz_file in self.recursive_files(indir):
            if gz_file.endswith('.inn.gz'):
                try:
                    res = ParseImgtAnnot(gz_file)()
                    data.append(res)
                    s += 1
                except Exception as e:
                    print(f"{gz_file}, error={e}")
                    f += 1
        outfile = os.path.join(self.data_dir, 'imgt_inn.json')
        self.save_json(data, outfile)
        print(f"{outfile}, succeed={s}, fail={f}")
    
    def build_inn_sequence(self, outprefix:str):
        res = {'chain':[],}
        seq_names = ['seq', 'seq_cdr', 'seq_cdr_mask', 'seq_align', \
            'seq_align_cdr','seq_align_cdr_mask',]
        for name in seq_names:
            key = f'vdomain_{name}'
            res[key] = []
    
        infile = os.path.join(self.data_dir, 'imgt_inn.json')
        data = self.load_data(infile, [])
        for rec in data:
            inn_number = rec['inn_number'][0]
            for chain in rec.get('chains', []):
                    chain_id = chain['chain_id']
                    if chain.get('chain_seq'):
                        desc = chain.get('chain_description')
                        res['chain'].append({
                            'inn': inn_number,
                            'id': chain_id,
                            'description': desc if desc else '',
                            'seq': chain['chain_seq'],
                        })
                    for domain in chain.get('v_domain', []):
                        domain_desc = domain.get('IMGT domain description', '')
                        for name in seq_names:
                            key = f"vdomain_{name}"
                            res[key].append({
                                'inn': inn_number,
                                'id': chain_id,
                                'description': domain_desc,
                                'seq': domain.get(name, '') ,
                            })

        for k,v in res.items():
            df = pd.DataFrame(v)
            df = df[['inn', 'id', 'description', 'seq']]
            outfile = f"{outprefix}_{k}.csv"
            df.to_csv(outfile, index=None)
            print(outfile)
        #save chain sequences to fasta
        fa_files = []
        for key in ('vdomain_seq',):
            outfile = f"{outprefix}_{key}.faa"
            ProcessSeq.to_fasta(res[key], outfile)
            fa_files.append(os.path.abspath(outfile))
        return fa_files


    '''
    download V-quest
    '''
    def download_vquest(self):
        info = {}
        endpoint = f"{self.url}download/V-QUEST/IMGT_V-QUEST_reference_directory/"
        response = requests.get(endpoint)
        species = re.findall(r'<a.*>(.+)/</a>', response.text)
        for specie in species:
            info[specie] = {}
            specie_url = f"{endpoint}{specie}/"
            response = requests.get(specie_url)
            specie_types = re.findall(r'<a.*>(.+)/</a>', response.text)
            for specie_type in specie_types:
                info[specie][specie_type] = []
                specie_type_url = f"{specie_url}{specie_type}/"
                response = requests.get(specie_type_url)
                file_names = re.findall(r'<a.*>(.+\.fasta)</a>', response.text)
                outdir = os.path.join(self.data_dir, 'V-QUEST', specie, specie_type)
                self.init_dir(outdir)
                for file_name in file_names:
                    file_url = os.path.join(specie_type_url, file_name)
                    outfile = os.path.join(outdir, file_name)
                    if not os.path.isfile(outfile):
                        with open(outfile, 'w') as f:
                            response = requests.get(file_url)
                            f.write(response.text)
                    info[specie][specie_type].append(outfile)
        return info

    '''
    download mAbDB
    '''
    def pull_mabdb(self, start:int=None, end:int=None):
        if start is None:
            start = 1 
        if end is None:
            end = 2000
        outdir = os.path.join(self.data_dir, 'mAb-DB')
        self.init_dir(outdir)

        s, f, k = 0 ,0, 0
        for abid in range(start, end):
            outfile = os.path.join(outdir, str(abid))
            if os.path.isfile(outfile):
                k += 1
            else:
                text = self.pull_mabdb_page(abid)
                if text and 'IMGT/mAb-DB ID' in text:
                    with open(outfile, 'w') as f:
                        f.write(text)
                    s += 1
                else:
                    f += 1
        print(f"Pull web page from IMGT/mAb-DB: succeed={s}, skipped={k}, failed={f}")

    def pull_mabdb_page(self, abid:str):
        try:
            endpoint = f"{self.url}/mAb-DB/mAbcard?AbId={abid}"
            response = requests.get(endpoint)
            if int(response.status_code) == 200:
                return response.text
        except Exception as e:
            pass
    '''
    retrieve mAbDB data
    '''
    def retrieve_mabdb(self, json_file):
        self.data = self.load_json(json_file)
        indir = os.path.join(self.data_dir, 'mAb-DB')
        text_iter = self.scan_text(indir)
        self.retrieve_data(text_iter, self.parse_imgt)
        self.retrieve_data(text_iter, self.parse_inn)
        self.retrieve_data(text_iter, self.parse_name)
        self.retrieve_data(text_iter, self.parse_product)
        self.retrieve_data(text_iter, self.parse_clinical)
        self.save_json(self.data, json_file)

    def parse_clinical(self, text) -> dict:
        res = {}
        pool = [
            ('clinical_trials', 'Clinical trials</td>'), 
            ('authority_decisions', 'Authority decisions</td>'),
            ('external_links', 'External links</td>'),
            ('imgt_notes', 'IMGT notes</td>'),
            ('biosimilars', 'Biosimilars</td>'),
        ]
        for key, tag in pool:
            substr = ParseXml.tag_substr(text, tag, '</td>')
            if  substr:
                res[key] = ParseXml.ul_li(substr)
        return res

    def parse_product(self, text) -> dict:
        res = {}
        pool = [
            ('expression_system', 'Expression system</td>'), 
            ('application', 'Application</td>'),
            ('clinical_domain', 'Clinical domain</td>'),
            ('mechanism_action', 'Mechanism of action</td>'),
            ('clinical_indication', 'Clinical indication</td>'),
            ('development_status', 'Development status</td>'),
        ]
        for key, tag in pool:
            substr = ParseXml.tag_substr(text, tag, '</td>')
            if  substr:
                dd = re.findall(r'<td.*>(.+)</td>', substr)
                if dd:
                    res[key] = dd

        substr = ParseXml.tag_substr(text, 'Company</td>', '</td>')
        if substr:
            res['company'] = ParseXml.href_text(substr)
        substr = ParseXml.tag_substr(text, 'Regulatory agency status and year</td>', '</td>')
        if substr:
            res['regulatory_status'] = ParseXml.ul_li(substr)
        return res

    def parse_imgt(self, text) -> dict:
        res = {}
        pool = [
            ('IMGT/mAb-DB ID</td>', 'imgt_mabdb'),
            ('IMGT/2Dstructure-DB</td>', 'imgt_2dstructure'),
            ('IMGT/3Dstructure-DB</td>', 'imgt_3dstructure'),
        ]
        for tag, key in pool:
            res[key] = {}
            substr = ParseXml.tag_substr(text, tag, '</td>')
            if  substr:
                _id = re.findall(r'<a.*>(.*)</a>', substr)
                if _id:
                    res[key]['id'] = _id[0]
                id_url = re.findall(r"href=\"(.*)\"", substr)
                if id_url:
                    res[key]['url'] =  f"{self.url}{id_url[0]}"
        
        res['engineered_variant'] = {}
        substr = ParseXml.tag_substr(text, 'IMGT engineered variant</td>', '</td>')
        if  substr:
            dd = re.findall(r'<td.*>(.+)</td>', substr)
            if dd:
                res['engineered_variant']['value'] = dd

        substr = ParseXml.tag_substr(text, 'IMGT engineered variant references</td>', '</td>')
        if  substr:
            res['engineered_variant']['references'] = ParseXml.ul_li(substr)
        return res

    def parse_name(self, text) -> dict:
        res = {}
        pool = [
            ('common_name', 'Common name</td>'), 
            ('proprietary_name', 'Proprietary name</td>'),
            ('species', 'Species</td>'),
            ('imgt_receptor_type', 'IMGT receptor type</td>'),
            ('receptor_identification', 'Receptor identification</td>'),
            ('conjugation','Radiolabelled / Conjugated / Fused</td>'),
        ]
        for key, tag in pool:
            substr = ParseXml.tag_substr(text, tag, '</td>')
            if  substr:
                dd = re.findall(r'<td.*>(.+)</td>', substr)
                if dd:
                    res[key] = dd
        return res

    def parse_inn(self, text) -> dict:
        key = 'inn'
        res = {key:{}}
        substr = ParseXml.tag_substr(text, 'INN</td>', '</td>')
        if  substr:
            dd = re.findall(r'<td.*>(.+)</td>', substr)
            if dd:
                res[key]['id'] = dd
        substr = ParseXml.tag_substr(text, 'INN Number</td>', '</td>')
        if  substr:
            dd = re.findall(r'<td.*>(.+)</td>', substr)
            if dd:
                res[key]['inn_number'] = dd
        substr = ParseXml.tag_substr(text, 'INN Prop. List</td>', '</td>')
        if  substr:
            dd = re.findall(r'<a.*>(.+)</a>', substr)
            if dd:
                res[key]['prop_list'] = dd
            dt = re.findall(r'\((\d+)\)', substr)
            if dt:
                res[key]['prop_list_year'] = dt
        substr = ParseXml.tag_substr(text, 'INN Rec. List</td>', '</td>')
        if  substr:
            dd = re.findall(r'<a.*>(.+)</a>', substr)
            if dd:
                res[key]['rec_list'] = dd
            dt = re.findall(r'\((\d+)\)', substr)
            if dt:
                res[key]['rec_list_year'] = dt
        return res

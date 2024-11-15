import sys
import os
import json

from bioomics import Uniref
from pipelines.utils import Utils

def main():
    n=0
    _m = ['UniProt', 'UniRef50', 'fasta']
    update_meta = {}
    meta_json = '/home/yuan/data/omics_data/epitope/index_meta.json'
    with open(meta_json, 'r') as f:
        meta = json.load(f)
    for k in list(meta):
        source = meta[k]['source']
        if _m not in source:
            update_meta[k] = meta[k]['json_file']

    print('parse sequence')
    n = 0
    uniprot_dir = '/home/yuan/data/UniProt/'
    uniref_iter = Uniref(uniprot_dir).scan_uniref50()
    keys = ['antigen', 'UniRef50']
    for rec in uniref_iter:
        if rec['id'] in update_meta:
            print(rec['id'])
            json_path = update_meta[rec['id']]
            Utils.update_json(json_path, rec, keys)
            meta[rec['id']]['source'].append(_m)
            if n%100 == 0:
                print(n, end=', ')
            n+=1
            # break
    with open(meta_json, 'w') as f:
        json.dump(meta, f, indent=4, sort_keys=True)
    # print(list(meta))
main()
import json
import os
import re
import requests
from biosequtils import Dir

class ProcessData:
    def __init__(self, data_dir:str):
        self.data_dir = data_dir
        self.data = None

    @staticmethod
    def init_dir(indir):
        '''
        create the directory if that doesn't exists
        '''
        pool = [indir,]
        while pool:
            curr_dir = pool.pop(0)
            if not os.path.isdir(curr_dir):
                parent_dir = os.path.dirname(curr_dir)
                if os.path.isdir(parent_dir):
                    try:
                        os.mkdir(curr_dir, 0o777)
                    except Exception as e:
                        print(e)
                        return False
                else:
                    pool = [parent_dir, curr_dir] + pool
        return True

    @staticmethod
    def recursive_files(indir): 
        '''
        list all files with a given directory and sub directories
        get all files
        '''
        for root, dirs, files in os.walk(indir):
            for filename in list(files):
                out_file = os.path.join(root, filename)
                if os.path.isfile(out_file) and out_file.find('/.') == -1:
                    yield out_file

    @staticmethod
    def load_json(infile):
        try:
            with open(infile, 'r') as f:
                return json.load(f)
        except Exception as e:
            pass
        return {}

    @staticmethod
    def save_json(data, outfile):
        try:
            with open(outfile, 'w') as f:
                json.dump(data, f, indent=4, sort_keys=True)
            return True
        except Exception as e:
            print(e)
        return False

    def load_data(self, json_file, default_data=None):
        # load data
        if os.path.isfile(json_file):
            self.data = ProcessData.load_json(json_file)
        else:
            self.data = default_data
        return self.data
    
    def retrieve_data(self, text_iter, func):
        n, s, k = 0, 0, 0
        for acc, text in text_iter:
            n += 1
            if acc not in self.data:
                self.data[acc] = {}
            # update
            val = func(text)
            if val:
                self.data[acc].update(val)
                s += 1
            else:
                k += 1
        print(f"Update data: {func}, succeed={s}, skipped={k}")

    def scan_text(self, indir:str=None):
        if indir is None:
            indir = self.data_dir
        file_iter = self.recursive_files(indir)
        for infile in file_iter:
            acc = os.path.basename(infile)
            with open(infile, 'r') as f:
                text = f.read()
                yield acc, text

    def pull_list(self, url):
        response = requests.get(url)
        names = re.findall(r'<a.*>(.+)</a>', response.text)
    

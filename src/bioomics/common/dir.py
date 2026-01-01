
import os
import re
import subprocess

class Dir:
    def __init__(self, indir):
        self.indir = indir

    def init_dir(self):
        '''
        create the directory if that doesn't exists
        '''
        pool = [self.indir,]
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

    def recursive_files(self): 
        '''
        list all files with a given directory and sub directories
        get all files
        '''
        for root, dirs, files in os.walk(self.indir):
            for filename in list(files):
                out_file = os.path.join(root, filename)
                if os.path.isfile(out_file) and out_file.find('/.') == -1:
                    yield out_file

    def compress_files(self, file_extension:str):
        '''
        compress files to *.gz
        '''
        info = {}
        for file_name in os.listdir(self.indir):
            path = os.path.join(self.indir, file_name)
            if path.endswith(file_extension):
                cmd = ['gzip', path]
                result = subprocess.run(cmd, capture_output=True, text=True)
                info[path].append(result)
        return info


    def format_dir(self):
        '''
        format directory
        '''
        #judge if absolute dir
        if self.indir.find('/')==0:
            pass
        elif self.indir.find('~')==0:
            self.indir = re.sub('~', os.getenv('HOME'), self.indir)
        elif self.indir is None:
            self.indir = os.getcwd()
        else: # ./test or test
            self.indir = os.path.abspath(self.indir)
        #alway followed by '/'
        if not self.indir[-1]=='/':
            self.indir += '/'
            
        #create directory if not exists
        if not os.path.isdir(self.indir):
            os.mkdir(self.indir, 0o777)
        return self.indir

    def clear_dir(self):
        '''
        delete all files
        '''
        all_files = self.recursive_files()
        for file in all_files:
            try:
                os.remove(file)
            except:
                pass

    def filter_recursive_files(self, extension:str):
        for path in self.recursive_files():
            if path.endswith(extension):
                yield path

    def recursive_dirs(self): 
        '''
        list all directories with a given directory and sub directories
        '''
        for root, dirs, files in os.walk(self.indir):
            if dirs == []:
                yield root

    @staticmethod
    def cascade_dir(parent_dir:str, id_str:str, num:int):
        '''
        restrict numbers of files or subdirectory for tuning up
        '''
        if len(id_str) < num:
            return parent_dir
        next_dir = os.path.join(parent_dir, id_str[:num])
        return Dir.cascade_dir(next_dir, id_str[num:], num)

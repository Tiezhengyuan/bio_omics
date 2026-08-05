'''
Execute third-party tools
default path ~/bio/tools/*
'''
from pathlib import Path
import subprocess

class Predictor:
    def __init__(self, args):
        self.args = args
    
    def _subprocess(self, cmd:list, shell:bool=False) -> bool:
        tool_name = self.args['tool_name']
        try:
            if self.args['verbose']:
                print(f"\n\nLoad {tool_name}: {' '.join(cmd)}")
            if shell:
                cmd = ' '.join(cmd)
            res = subprocess.run(cmd, shell=shell, capture_output=True, text=True)
            if res.returncode == 0:
                return True
            if self.args['verbose']:
                print('status: ', res.returncode)
                print('stdout: ', res.stdout)
                print('stderr: ', res.stderr)
        except Exception as e:
            print(f"ERROR: failed in execution of {tool_name}, error={e}")
        return False
    
    def colabfold_batch(self) -> bool:
        '''
        AlphaFold2: predict protein structure given sequences
        '''
        if 'tool_name' not in self.args:
            self.args['tool_name'] = 'ColabFold_Batch'
        tool_dir = self.args['tools_dir'] / 'localcolabfold'
        bin_dir = tool_dir / Path('.pixi/envs/default/bin')
        
        # I. MSA locally
        # cmd0 = [
        #     'pixi', 'run',
        #     '--manifest-path', str(tool_dir / 'pyproject.toml'),
        #     str(bin_dir / 'colabfold_batch'),
        #     '--msa-only',
        #     str(self.args['faa_file']),
        #     str(self.args['outdir']),
        # ]
        # self._subprocess(cmd0)
                
        # II. predict
        cmd = [
            'pixi', 'run', 
            '--manifest-path', str(tool_dir / 'pyproject.toml'),
            str(bin_dir / 'colabfold_batch'),
            '--templates', '--random-seed', '42',
        ]
        if self.args.get('save_all', False):
            cmd.append('--save-all')
        if self.args.get('relax', True):
            cmd += ['--amber', '--use-gpu-relax',]
            if 'max_iterations' in self.args:
                cmd += [
                    '--relax-max-iterations', str(self.args['max_iterations']),
                ]
        cmd += [
            '--model-type', str(self.args.get('model_type', 'auto')),
            '--num-recycle', str(self.args.get('num_recycle', 3)),
            str(self.args['faa_file']),
            str(self.args['outdir']),
        ]
        return self._subprocess(cmd)

    def haddock3(self):
        '''
        haddock3: docking simulation
        '''
        if 'tool_name' not in self.args:
            self.args['tool_name'] = 'HadDock3'
        cmd = ['haddock3', str(self.args['cfg_file']),]
        return self._subprocess(cmd)

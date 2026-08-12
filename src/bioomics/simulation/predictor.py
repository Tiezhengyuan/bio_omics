'''
Execute third-party tools
default path ~/bio/tools/*
'''
import asyncio
from pathlib import Path
import subprocess

class Predictor:
    def __init__(self, input, outdir:Path, tools_dir:Path=None, verbose:bool=False):
        self.input = input
        self.outdir = outdir
        self.outdir.mkdir(parents=True, exist_ok=True)
        self.tools_dir = tools_dir
        self.tool_name = None
        self.verbose = verbose
    
    @staticmethod
    async def run_command(cmd, tool_name):
        print(f"### {tool_name} asynchronous running: {cmd}")
        process = await asyncio.create_subprocess_exec(
            *cmd,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
        )
        stdout, stderr = await process.communicate()
        returncode = process.returncode
        if returncode != 0:
            print(stderr.decode())
            return False
        print(stdout.decode())
        return True

    @staticmethod
    def run_subprocess(cmd:list, tool_name:str=None, shell:bool=False) -> bool:
        try:
            print(f"\n\nLoad {tool_name}: {' '.join(cmd)}")
            if shell:
                cmd = ' '.join(cmd)
            res = subprocess.run(cmd, shell=shell, capture_output=True, text=True)
            print('status: ', res.returncode)
            print('stdout: ', res.stdout)
            print('stderr: ', res.stderr)
            return res
        except Exception as e:
            print(f"ERROR: failed in {tool_name} execution, error={e}")
        return None
    
    
    @staticmethod
    def colabfold_search(params:dict):
        '''
        AlphaFold2: local MSA
        params = {
            tools_dir:Path,
            query_faa:Path,
            outdir:Path,
            colabfold_db:Path
        }
        '''
        tool_name = 'colabfold_search'
        cmd = [
            'conda', 'run', 
            '-n', 'colabfold',
            '--cwd', str(params['tools_dir'] / 'ColabFold'),
            tool_name,
            '--gpu', '1',
            '--threads', '8',
            '--use-env', '0',
            '--use-templates', '1',
            '--db1', 'uniref30_2302_db',
            '--db2', 'pdb100_230517',
            '--db3', 'colabfold_envdb_202108_db',
            str(params['query_faa']),
            str(params['colabfold_db']),
            str(params['outdir']),
        ]
        return Predictor.run_subprocess(cmd, tool_name)

    @staticmethod
    def colabfold_batch_msa(params:dict):
        '''
        AlphaFold2: access public MSA server
        params = {
            tools_dir: Path,
            query_faa: Path,
            outdir: Path,
        }
        '''
        tool_name = 'colabfold_batch'
        cmd = [
            'conda', 'run', '-n', 'colabfold',
            '--cwd', str(params['tools_dir'] / 'ColabFold'),
            tool_name,
            '--msa-only',
            '--templates',
            str(params['query_faa']),
            str(params['outdir']),
        ]
        return Predictor.run_subprocess(cmd, tool_name)

    @staticmethod
    async def colabfold_batch(params:dict):
        '''
        AlphaFold2: predict protein structure given sequences
        params = {
            tools_dir: Path,
            input: Path,
            outdir: Path,
            ...
        }
        '''
        tool_name = 'colabfold_batch'
        cmd = [
            'conda', 'run',
            '-n', 'colabfold',
            '--cwd', str(params['tools_dir'] / 'ColabFold'),
            tool_name,
            '--num-models', '5', \
            '--model-order', '1,2,3,4,5', \
            '--random-seed', '42',
        ]
        if params.get('use_templates'):
            cmd.append('--templates')
        if params.get('save_all'):
            cmd.append('--save-all')
        if params.get('relax'):
            cmd += ['--amber', '--use-gpu-relax',]
            if 'max_iterations' in params:
                cmd += [
                    '--relax-max-iterations', str(params['max_iterations']),
                ]
        cmd += [
            '--model-type', str(params.get('model_type', 'auto')),
            '--num-recycle', str(params.get('num_recycle', 3)),
            str(params['input']), #fasta or a3m
            str(params['outdir']),
        ]
        return Predictor.run_subprocess(cmd, tool_name)


################################################################################3

    @staticmethod
    def haddock3(params:dict):
        '''
        haddock3: docking simulation
        '''
        tool_name = 'HadDock3'
        cmd = ['haddock3', str(params['cgi_file']),]
        return Predictor._subprocess(cmd, False, tool_name)

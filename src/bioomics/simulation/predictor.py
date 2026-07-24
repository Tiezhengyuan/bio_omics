import os
import subprocess

class Predictor:
    def __init__(self, args):
        self.args = args
    
    def colabfold(self):
        cmd = ['colabfold_batch', '--templates', '--random-seed', '42',]
        if self.args.get('save_all', False):
            cmd.append('--save-all')
        if self.args.get('relax', True):
            cmd += ['--amber', '--use-gpu-relax',]
            if 'max_iterations' in self.args:
                cmd += ['--relax-max-iterations', self.args['max_iterations'],]
        cmd += [
            '--model-type', self.args.get('model_type', 'auto'),
            '--num-recycle', self.args.get('num_recycle'),
            self.args['faa_file'], self.args['output_dir'],
        ]
        try:
            if self.args['verbose']:
                print('\n\nLoad colabfold: ', ' '.join(cmd))
            # cmd = ' '.join(cmd)
            result = subprocess.run(cmd, shell=False,
                capture_output=True, text=True)
            if self.args['verbose']:
                print('status: ', result.returncode)
                print('stdout: ', result.stdout)
                print('stderr: ', result.stderr)
            return True
        except Exception as e:
            print('ERROR: failed in alphafold through colabfold, error=', e)
        return False

    def haddock3(self):
        cmd = ['haddock3', self.args['cfg_file'],]
        try:
            if self.args['verbose']:
                print('\n\nLoad haddock3: ', ' '.join(cmd))
            # cmd = ' '.join(cmd)
            result = subprocess.run(cmd, shell=False,
                capture_output=True, text=True)
            if self.args['verbose']:
                print('status: ', result.returncode)
                print('stdout: ', result.stdout)
                print('stderr: ', result.stderr)
            return True
        except Exception as e:
            print('ERROR: failed in haddock3, error=', e)
        return False

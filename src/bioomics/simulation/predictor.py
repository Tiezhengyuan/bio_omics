import os
import subprocess

class Predictor:
    def __init__(self, args):
        self.args = args
    
    def colabfold(self):
        cmd = ['colabfold_batch', 
            '--templates', '--save-all',
            '--amber', '--use-gpu-relax',
            '--random-seed', '42',
            '--model-type', self.args.get('model_type', 'auto'),
            '--relax-max-iterations', self.args.get('max_iterations', '2000'), 
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

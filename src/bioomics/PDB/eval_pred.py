import os
import pandas as pd
from pyrosetta import pyrosetta, io, rosetta
import subprocess

class EvalPred:

    def __init__(self, bin_dir=None, data_dir=None):
        self.bin_dir = bin_dir
        self.data_dir = data_dir

    def rmsd(self, native_pdb, model_pdb):
        try:
            pyrosetta.init()
            native_pose = io.pose_from_pdb(native_pdb)
            model_pose = io.pose_from_pdb(model_pdb)        
            rmsd = rosetta.core.scoring.CA_rmsd(
                native_pose, model_pose
            )
            return rmsd
        except Exception as e:
            pass
        return None

    def TMalign(self, pdb1, pdb2):
        try:
            exe = os.path.join(self.bin_dir, 'MMalign')
            cmd = [exe, pdb1, pdb2, '-outfmt', '2',]
            result = subprocess.run(cmd, capture_output=True,
                text=True, check=True)
            # print("STDOUT:", result.stdout)
            # print("STDERR:", result.stderr)
            lines = result.stdout.split('\n')
            l0 = lines[0].split('\t')
            l1 = lines[1].split('\t')
            res = {k:v for k,v in zip(l0, l1)}
            return res
        except Exception as e:
            print(e)
            pass
        return {}
        
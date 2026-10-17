import os
from pathlib import Path
from pyrosetta import pyrosetta, io, rosetta
import subprocess

class EvalPred:

    def __init__(self, bin_dir:Path=None):
        self.bin_dir = bin_dir

    def rmsd(self, native_pdb:str, model_pdb:str):
        try:
            pyrosetta.init()
            native_pose = io.pose_from_pdb(str(native_pdb))
            model_pose = io.pose_from_pdb(str(model_pdb))
            rmsd = rosetta.core.scoring.CA_rmsd(
                native_pose, model_pose
            )
            return rmsd
        except Exception as e:
            pass
        return None

    def batch_rmsd(self, pool:list):
        pyrosetta.init()
        for row in pool:
            try:
                native_pose = io.pose_from_pdb(str(row['pdb1']))
                model_pose = io.pose_from_pdb(str(row['pdb2']))
                row['rmsd'] = rosetta.core.scoring.CA_rmsd(
                    native_pose, model_pose
                )
                yield row
            except Exception as e:
                pass

    
    def TMalign(self, pdb1:str, pdb2:str):
        try:
            exe = str(self.bin_dir / 'MMalign')
            cmd = [exe, str(pdb1), str(pdb2), '-outfmt', '2',]
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
    
    def prodigy(self, pdb_file:Path, outfile:Path=None):
        try:
            cmd = ['prodigy', str(pdb_file),]
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            if outfile is not None:
                with open(outfile, 'w') as f:
                    f.write(result.stdout)
            # parse result
            res = self.parse_prodigy(result.stdout)
            return res
        except Exception as e:
            print('error=', e)
        return None

    def parse_prodigy(self, text):
        # not contacts
        if 'No contacts found' in text:
            return {'contacts': False, }
        
        res = {'contacts': True, }
        for line in text.split('\n'):
            if 'No.' in line or 'Percentage' in line or 'Predicted' in line:
                k, v = line.split(':')
                val = float(v.strip())
                if 'intermolecular' in line:
                    res['intermolecular_contacts'] = val
                elif 'charged-charged' in line:
                    res['charged-charged_contacts'] = val
                elif 'charged-polar' in line:
                    res['charged-polar_contacts'] = val
                elif 'charged-apolar' in line:
                    res['charged-apolar_contacts'] = val
                elif 'polar-polar' in line:
                    res['polar-polar_contacts'] = val
                elif 'apolar-polar' in line:
                    res['apolar-polar_contacts'] = val
                elif 'apolar-apolar' in line:
                    res['apolar-apolar_contacts'] = val
                elif 'apolar NIS' in line:
                    res['apolar_NIS_residues'] = val
                elif 'charged NIS' in line:
                    res['charged_NIS_residues'] = val
                elif 'binding affinity' in line:
                    res['binding_affinity'] = val
                elif 'dissociation constant' in line:
                    res['dissociation_constant'] = val
        return res
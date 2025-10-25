from pyrosetta import pyrosetta, io, rosetta

class AnalyzeModel:

    @staticmethod
    def rmsd(native_pdb, model_pdb):
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


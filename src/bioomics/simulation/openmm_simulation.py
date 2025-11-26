'''
MD simulation using openMM
'''
import gzip
from openmm import app, openmm, unit
import os
from sys import stdout
from ..common import ProcessJson

FORCE_FIELDS = {
    # a comprehensive set of Amber14 force field files
    'general': ['amber14-all.xml',],
    'protein': ['amber14/protein.ff14SB.xml', 'amber99sb.xml',],
    'membrane_protein': ['charmm36.xml',],
    'nucleic': ['amber14/all_nucleic.xml',],
    'DNA': ['amber14/DNA.OL15.xml',],
    'RNA': ['amber14/RNA.OL3.xml',],
    # explicit solvent: water
    'explicit_solvent': [
        'amber14/tip3p.xml',    # most common, standard TIP3P
        'amber14/tip3pfb.xml',  # updated version of tip3p
        'amber14/spce.xml',     # SPC/E
        'amber14/opc.xml',      # higher accuracy
    ],
    'implicit_solvent': ['implicit/gbn2.xml', 'implicit/obc.xml',],
    'small_molecule': ['gaff.xml', 'gaff2.xml', 'charmm36/chem.xml',],
    # combined with explicit solvent
    'ion': ['amber14/ions.xml',],
}

# water system
COMPLEX_FF = {
    'general_use': ['amber14-all.xml', 'amber14/tip3p.xml',],
    'protein_protein': ['amber14/protein.ff14SB.xml', 'amber14/tip3pfb.xml',],
    'membrane_protein': ['charmm36.xml', 'amber14/tip3pfb.xml',],
    'channel_protein': ['charmm36.xml', 'amber14/ions.xml', 'amber14/tip3pfb.xml',],
    'protein_dna': ['amber14/protein.ff14SB.xml', 'amber14/DNA.OL15.xml'],
    'protein_rna': ['amber14/protein.ff14SB.xml', 'amber14/RNA.OL13.xml'],
    'protein_ion': ['amber14/protein.ff14SB.xml', 'amber14/ions.xml', 'amber14/tip3pfb.xml',],
    'protein_ligand': ['amber14/protein.ff14SB.xml', 'gaff.xml', 'amber14/tip3pfb.xml',],
}

class OpenmmSimulation:

    def __init__(self, pdb_file, outdir:str=None, verbose:bool=False):
        self.pdb_file = pdb_file
        self.outdir = outdir
        self.verbose = verbose
        self.pdb = None
        self.simulation = None 

    def __call__(self):
        self.info = {}
        self.load_pdb()
        if self.pdb:
            self.create_system()
            # energy minimization
            num_state = 10
            max_iterations=0
            report_interval = 1000
            self.minimize_energy(num_state, max_iterations, report_interval)
            return self.export_state()
        return 'fail'

    def load_pdb(self):
        '''
        create self.pdb
        '''
        try:
            if self.pdb_file.endswith('.gz'):
                with gzip.open(self.pdb_file, 'rt') as f:
                    self.pdb = app.PDBFile(pdbfile=f)
            else:
                self.pdb = app.PDBFile(self.pdb_file)
        except Exception as e:
            print(f"Error: wrong pdb file: {self.pdb_file}")
            self.info['error_load_pdb'] = str(e)

    def create_system(self):
        '''
        create self.simulation
        '''
        ff_model = 'amber19-all.xml'
        water_model = 'amber19/tip3pfb.xml'
        try:
            forcefield = app.ForceField(ff_model, water_model)
            system = forcefield.createSystem(
                self.pdb.topology,
                nonbondedMethod = app.PME, # PME periodic, Nocutoff
                nonbondedCutoff = 1 * unit.nanometer,
                constraints = app.HBonds, # HBonds, AllBonds, None
                rigidWater=True        # Keep water molecules rigid
            )
            self.info['force_field'] = [ff_model, water_model]
            # Langevin dynamics, with the LFMiddle discretization
            # (J. Phys. Chem. A 2019, 123, 28, 6056-6079).
            integrator = openmm.LangevinMiddleIntegrator(
                300 * unit.kelvin, # 300K is room temperature
                1.0 / unit.picosecond,
                0.002 * unit.picoseconds, # 0.001-0.004 ps for stability
            )

            self.simulation = app.Simulation(self.pdb.topology, system, integrator)
            self.simulation.context.setPositions(self.pdb.positions)
            curr_state = self.simulation.context.getState(getEnergy=True)
            self.info['original_energy'] = str(curr_state.getPotentialEnergy())
            if self.verbose:
                print("Initial energy:", self.info['original_energy'])
        except Exception as e:
            print(f"Warning: System creation failed {self.pdb_file}")
            self.info['error_create_system'] = str(e)

    def minimize_energy(self, num_state:int, max_iterations:int, report_interval:int):
        '''
        create self.simulation
        '''
        try:
            # energy minimization
            self.simulation.minimizeEnergy(maxIterations=max_iterations)

            if os.path.isdir(self.outdir):
                # checkpoint
                outfile = os.path.join(self.outdir, 'checkpoint.chk')
                self.simulation.reporters.append(
                    app.CheckpointReporter(outfile, report_interval)
                )
                # dynamic coordinate data
                outfile = os.path.join(self.outdir, 'trajactory.dcd')
                self.simulation.reporters.append(
                    app.DCDReporter(outfile, report_interval)
                )
                # pdb at each state
                outfile = os.path.join(self.outdir, 'structure.pdb')
                self.simulation.reporters.append(
                    app.PDBReporter(outfile, report_interval)
                )
                # meta data of simulation
                outfile = os.path.join(self.outdir, 'state_data.csv')
                self.simulation.reporters.append(
                    app.StateDataReporter(outfile, report_interval, step=True, \
                        potentialEnergy=True, temperature=True)
                )
            else:
                self.simulation.reporters.append(
                    app.StateDataReporter(stdout, report_interval, step=True, \
                        potentialEnergy=True, temperature=True)
                )

            steps = report_interval * num_state
            self.simulation.step(steps)
            minmized_state = self.simulation.context.getState(getEnergy=True)
            self.info["after_minimization"] = str(minmized_state.getPotentialEnergy())
            if self.verbose:
                print("Energy after minimization:", self.info["after_minimization"])
        except Exception as e:
            print(f"Warning: energy minimization failed {self.pdb_file}")
            self.info['error_energy_minimization'] = str(e)

    def export_state(self):
        if not os.path.isdir(self.outdir):
            if self.verbose:
                print('No exports', self.info)
        print(self.info)
        try:       
            # Save the current state of the simulation to a file.
            outfile = os.path.join(self.outdir, 'relaxed_state.xml')
            self.simulation.saveState(outfile)

            # Save minimized structure
            relaxed_pdb = os.path.join(self.outdir, 'relaxed.pdb')
            self.info['relaxed_pdb'] = relaxed_pdb
            positions = self.simulation.context.getState(getPositions=True).getPositions()
            with open(relaxed_pdb, 'w') as f:
                app.PDBFile.writeFile(self.simulation.topology, positions, f)
            
            # save self.info
            json_file = os.path.join(self.outdir, 'info.json')
            ProcessJson.save_json(self.info, json_file)
            return 'succeed'
        except Exception as e:
            print(f"Warning: No export energy minimization. {self.pdb_file}")
            self.info['error_save_state'] = str(e)
        return 'fail'

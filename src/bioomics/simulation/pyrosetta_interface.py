'''
scoring of designed binders using Rosetta interface metrics for running PyRosetta within BindCraft
that would be called within the main design loop in run_design_loop.py 

# backward-compatible procedural wrapper
def score_interface(pdb_file, binder_chain="B"):
    return PyRosettaInterface().run(pdb_file, binder_chain)

prompt for unit testing:
    create or add unit testing cases 
        - sourcde code script is denovo/filter/pyrosetta_interface.py
        - DONOT change source code unless there is ANY bug.
        - testing cases are saved into tests/test_pyrosetta_interface.py.
        - input data are put in tests/data.

prompt for legacy coding updates:
    analyze denovo/filterpyrosetta_interface.py
    - debugging if there is
    - optimize for speeding up execution.
'''
import pyrosetta as pr
from pyrosetta.rosetta.core.select.residue_selector import ChainSelector
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects

try:
    # available in pyrosetta versions that expose the DockingPartners API
    from pyrosetta.rosetta.core.pose import DockingPartners
except ImportError:
    DockingPartners = None

from denovo.filter import Trajectory


class PyRosettaInterface:
    def __init__(self, dalphaball_path:str, scorefxn=None):
        pr.init(f"-mute all -holes:dalphaball {dalphaball_path}")
        self.scorefxn = scorefxn or pr.get_fa_scorefxn()
                    
    def collect_interface_residues(self, pdb_file, binder_chain):
        '''return (interface_AA, interface_residues) of the binder interface'''
        self.interface_scores = {}
        interface_residues_set = Trajectory(pdb_file).hotspot_residues(binder_chain)
        interface_AA = {aa: 0 for aa in "ACDEFGHIKLMNPQRSTVWY"}
        interface_residues = []
        for pdb_res_num, aa_type in interface_residues_set.items():
            interface_AA[aa_type] += 1
            interface_residues.append(f"{binder_chain}{pdb_res_num}")

        # interface_nres: Number of residues at the interface, defined as residues with any atom within 5 Å
        interface_nres = len(interface_residues)
        self.interface_scores["interface_nres"] = interface_nres

        return interface_AA, interface_residues
    
    def compute_interface_scores(self, pdb_file:str, binder_chain="B",
            interface_AA=None, interface_residues=None) -> dict:
        '''
        score pdb_file and return the rounded interface_scores dict;
        interface_AA / interface_residues default to a fresh
        collect_interface_residues call when not provided
        '''
        if interface_AA is None or interface_residues is None:
            interface_AA, interface_residues = self.collect_interface_residues(pdb_file, binder_chain)

        pose = pr.pose_from_pdb(pdb_file)
        iam = self._run_interface_analyzer(pose)
        interfacescore = iam.get_all_data()
        
        # interface_hydrophobicity
        self._interface_hydrophobicity(interface_AA)
        # binder_sore
        self._binder_score(pose, binder_chain)
        # surface_hydrophobicity
        self._surface_hydrophobicity(pose, binder_chain)
        # binder_sasa
        binder_sasa = self._binder_sasa(pose, binder_chain)
        # interface_dSASA
        interface_dSASA = self._interface_dSASA(iam)
        # interface_fraction
        self._interface_binder_fraction(binder_sasa, interface_dSASA)
        # interface_dG
        self._interface_dG(iam)
        #interface_packstat
        self._interface_packstat(iam)
        # interface_sc
        self._interface_sc(interfacescore)
        # interface_interface_hbonds
        interface_interface_hbonds = self._interface_interface_hbonds(interfacescore)
        # interface_dG_SASA_ratio
        self._interface_dG_SASA_ratio(interfacescore)
        # interface_hbond_percentage
        self._interface_hbond_percentage(interface_interface_hbonds)
        # interface_delta_unsat_hbonds
        interface_delta_unsat_hbonds = self._interface_delta_unsat_hbonds(pose)
        # interface_deleta_unsat_hbonds_percentage
        self._interface_bunsch_percentage(interface_delta_unsat_hbonds)

        return {k: round(v, 2) if isinstance(v, float) else v
            for k, v in self.interface_scores.items()
        }


    def _run_interface_analyzer(self, pose):
        # Run Rosetta InterfaceAnalyzerMover on the given pose to compute interface metrics. 
        iam = InterfaceAnalyzerMover()
        if DockingPartners is not None:
            partners = DockingPartners.docking_partners_from_string("A_B")
            iam.set_interface(partners)
        else:
            # used by old version of pyrosetta
            iam.set_interface("A_B")
        iam.set_scorefunction(self.scorefxn)
        iam.set_compute_packstat(True)
        iam.set_compute_interface_energy(True)
        iam.set_calc_dSASA(True)
        iam.set_calc_hbond_sasaE(True)
        iam.set_compute_interface_sc(True)
        iam.set_pack_separated(True)
        iam.apply(pose)
        return iam

    def _interface_hydrophobicity(self, interface_AA):
        interface_nres = self.interface_scores['interface_nres']
        # Percentage of hydrophobic residues at the interface based on the amino acid composition. 
        hydrophobic_aa = set("ACFILMPVWY")
        hydrophobic_count = sum(interface_AA[aa] for aa in hydrophobic_aa)
        interface_hydrophobicity = (hydrophobic_count / interface_nres) * 100 if interface_nres else 0
        self.interface_scores["interface_hydrophobicity"] = interface_hydrophobicity
        return interface_hydrophobicity

    def _interface_dSASA(self, iam):
        # Interface delta SASA (dSASA)
        interface_dSASA = iam.get_interface_delta_sasa()
        self.interface_scores["interface_dSASA"] = interface_dSASA
        return interface_dSASA

    def _interface_dG(self, iam):
        # Interface energy (dG)
        interface_dG = iam.get_interface_dG()
        self.interface_scores["interface_dG"] = interface_dG
        return interface_dG
    
    def _interface_packstat(self, iam):
        # Interface packing statistic (packstat)
        interface_packstat = iam.get_interface_packstat()
        self.interface_scores["interface_packstat"] = interface_packstat
        return interface_packstat
    
    def _interface_sc(self, interfacescore):
        interface_sc = interfacescore.sc_value
        self.interface_scores["interface_sc"] = interface_sc
        return interface_sc
    
    def _interface_interface_hbonds(self, interfacescore):
        interface_interface_hbonds = interfacescore.interface_hbonds
        self.interface_scores["interface_interface_hbonds"] = interface_interface_hbonds
        return interface_interface_hbonds
    
    def _interface_dG_SASA_ratio(self, interfacescore):
        interface_dG_SASA_ratio = interfacescore.dG_dSASA_ratio * 100
        self.interface_scores["interface_dG_SASA_ratio"] = interface_dG_SASA_ratio
        return interface_dG_SASA_ratio

    def _interface_hbond_percentage(self, interface_interface_hbonds):
        interface_nres = self.interface_scores['interface_nres']
        # Percentage of interface residues that are involved in hydrogen bonds at the interface.
        interface_hbond_percentage = (interface_interface_hbonds / interface_nres) * 100 if interface_nres else None
        self.interface_scores["interface_hbond_percentage"] = interface_hbond_percentage
        return interface_hbond_percentage

    def _interface_delta_unsat_hbonds(self, pose):
        buns_filter = XmlObjects.static_get_filter(
            '<BuriedUnsatHbonds report_all_heavy_atom_unsats="true" scorefxn="scorefxn" ignore_surface_res="false" use_ddG_style="true" dalphaball_sasa="1" probe_radius="1.1" burial_cutoff_apo="0.2" confidence="0" />'
        )
        interface_delta_unsat_hbonds = buns_filter.report_sm(pose)
        self.interface_scores["interface_delta_unsat_hbonds"] = interface_delta_unsat_hbonds
        return interface_delta_unsat_hbonds

    def _interface_bunsch_percentage(self, interface_delta_unsat_hbonds):
        interface_nres = self.interface_scores['interface_nres']
         # Percentage of interface residues that are buried unsatisfied hydrogen bonds at the interface.
        interface_bunsch_percentage = (interface_delta_unsat_hbonds / interface_nres) * 100 if interface_nres else None
        self.interface_scores["interface_delta_unsat_hbonds_percentage"] = interface_bunsch_percentage
        return interface_bunsch_percentage

    def _binder_score(self, pose, binder_chain):
        # Total score of the binder chain            
        chain_design = ChainSelector(binder_chain)
        tem = pr.rosetta.core.simple_metrics.metrics.TotalEnergyMetric()
        tem.set_scorefunction(self.scorefxn)
        tem.set_residue_selector(chain_design)
        binder_score = tem.calculate(pose)
        self.interface_scores["binder_score"] = binder_score
        return binder_score

    def _binder_sasa(self, pose, binder_chain):
        # SASA of the binder chain
        chain_design = ChainSelector(binder_chain)
        bsasa = pr.rosetta.core.simple_metrics.metrics.SasaMetric()
        bsasa.set_residue_selector(chain_design)
        binder_sasa = bsasa.calculate(pose)
        self.interface_scores["binder_sasa"] = binder_sasa
        return binder_sasa
    
    #interface_binder_fraction renamed to interface_fration in interface_scores
    def _interface_binder_fraction(self, binder_sasa, interface_dSASA):
        # Fraction of the binder surface that is at the interface (interface_dSASA / binder_sasa)
        interface_binder_fraction = (interface_dSASA / binder_sasa) * 100 if binder_sasa > 0 else 0
        self.interface_scores["interface_fraction"] = interface_binder_fraction
        return interface_binder_fraction

    def _surface_hydrophobicity(self, pose, binder_chain):
        binder_pose = {
            pose.pdb_info().chain(pose.conformation().chain_begin(i)): p
            for i, p in zip(range(1, pose.num_chains() + 1), pose.split_by_chain())
        }[binder_chain]

        layer_sel = pr.rosetta.core.select.residue_selector.LayerSelector()
        layer_sel.set_layers(pick_core=False, pick_boundary=False, pick_surface=True)
        surface_res = layer_sel.apply(binder_pose)

        exp_apol_count = 0
        total_count = 0
        for i in range(1, len(surface_res) + 1):
            if surface_res[i]:
                res = binder_pose.residue(i)
                if res.is_apolar() or res.name() in ("PHE", "TRP", "TYR"):
                    exp_apol_count += 1
                total_count += 1

        surface_hydrophobicity = exp_apol_count / total_count if total_count else 0
        self.interface_scores["surface_hydrophobicity"] = surface_hydrophobicity
        return surface_hydrophobicity

    def run(self, pdb_file:str, binder_chain="B") -> tuple:
        # interface_AA, interface_residues
        interface_AA, interface_residues = self.collect_interface_residues(pdb_file, binder_chain)
        # interface_scores
        interface_scores = self.compute_interface_scores(
            pdb_file, binder_chain, interface_AA, interface_residues)
        return interface_AA, interface_residues, interface_scores

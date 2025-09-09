"""
Fix and minimize a protein structure.
"""
import os
import tempfile
from io import StringIO
from pdbfixer import PDBFixer
from openmm import app, unit
import openmm as mm
from Bio.PDB import PDBIO

def biopy_structure_to_pdb_text(structure):
    """
    Convert a Biopython structure to a PDB file string.
    """
    io = PDBIO()
    io.set_structure(structure)
    sio = StringIO()
    io.save(sio)
    return sio.getvalue()

class FixerMinimizer:
    """
    A class to fix and minimize a protein structure.
    """
    def __init__(self, ph=7.0, platform_name=None):
        """
        Initialize the FixerMinimizer.
        """
        self.ph = ph
        self.platform_name = platform_name

    def fix_and_minimize(self, biopy_struct, keep_water=True, niter_minimize=500, restrain_calpha=True):
        """
        Fix and minimize a protein structure.
        """
        pdb_text = biopy_structure_to_pdb_text(biopy_struct)
        sio = StringIO(pdb_text)
        fixer = PDBFixer(pdbfile=sio)
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(self.ph)
        if not keep_water:
            fixer.removeHeterogens(keepWater=False)
        topology = fixer.topology
        positions = fixer.positions

        forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3p.xml')
        if topology.getPeriodicBoxVectors():
            system = forcefield.createSystem(topology,
                                             nonbondedMethod=app.PME,
                                             nonbondedCutoff=1.0*unit.nanometer,
                                             constraints=app.HBonds)
        else:
            system = forcefield.createSystem(topology,
                                             nonbondedMethod=app.NoCutoff,
                                             constraints=app.HBonds)

        if restrain_calpha:
            from openmm import CustomExternalForce
            k_rest = 10.0 * unit.kilocalories_per_mole / (unit.angstrom**2)
            restraint = CustomExternalForce("0.5 * k * ((x - x0)^2 + (y - y0)^2 + (z - z0)^2)")
            restraint.addPerParticleParameter("x0")
            restraint.addPerParticleParameter("y0")
            restraint.addPerParticleParameter("z0")
            restraint.addGlobalParameter("k", k_rest)
            for atom in topology.atoms():
                if atom.name == "CA":
                    idx = atom.index
                    pos = positions[idx]
                    restraint.addParticle(idx, [pos.x, pos.y, pos.z])
            system.addForce(restraint)

        integrator = mm.LangevinIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 0.002*unit.picoseconds)
        platform = None
        if self.platform_name:
            platform = mm.Platform.getPlatformByName(self.platform_name)
        simulation = app.Simulation(topology, system, integrator, platform) if platform else app.Simulation(topology, system, integrator)
        simulation.context.setPositions(positions)

        simulation.minimizeEnergy(maxIterations=niter_minimize)
        state = simulation.context.getState(getPositions=True, enforcePeriodicBox=False)
        minimized_positions = state.getPositions()

        with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", mode="w", encoding="utf-8") as out_pdb_path:
            app.PDBFile.writeFile(simulation.topology, minimized_positions, out_pdb_path)
            temp_name = out_pdb_path.name

        with open(temp_name, 'r', encoding="utf-8") as f:
            out_pdb_text = f.read()
        os.remove(temp_name)
        return out_pdb_text
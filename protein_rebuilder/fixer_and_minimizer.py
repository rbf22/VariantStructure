# protein_rebuilder/fixer_and_minimizer.py
from pdbfixer import PDBFixer
from openmm import app, unit, LangevinIntegrator
import openmm as mm
import tempfile
import os
from Bio.PDB import PDBIO

def biopy_structure_to_pdb_text(structure):
    """Return PDB text for a Bio.PDB structure object (string)."""
    fh = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb")
    tmpname = fh.name
    fh.close()
    io = PDBIO()
    io.set_structure(structure)
    io.save(tmpname)
    with open(tmpname, 'r') as f:
        pdb_text = f.read()
    os.remove(tmpname)
    return pdb_text

class FixerMinimizer:
    """
    Uses PDBFixer to add missing atoms & sidechains, add hydrogens, then
    minimize using OpenMM.
    """

    def __init__(self, ph=7.0, platform_name=None):
        self.ph = ph
        self.platform_name = platform_name

    def fix_and_minimize(self, biopy_struct, keep_water=True, niter_minimize=500):
        pdb_text = biopy_structure_to_pdb_text(biopy_struct)
        # create a tempfile with the pdb for PDBFixer
        fixer = PDBFixer(pdbfile=None)
        # PDBFixer can be instantiated empty then readPDB...
        fixer.readPDB(pdb_text)
        # Find missing
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        # optionally add missing hydrogens
        fixer.addMissingHydrogens(self.ph)
        # Optionally remove or keep water
        if not keep_water:
            fixer.removeHeterogens(keepWater=False)
        # Now create OpenMM objects
        # Convert to OpenMM topology+positions
        topology = fixer.topology
        positions = fixer.positions

        forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3p.xml')
        system = forcefield.createSystem(topology,
                                         nonbondedMethod=app.PME,
                                         nonbondedCutoff=1.0*unit.nanometer,
                                         constraints=app.HBonds)

        # optional: add a positional restraint on backbone heavy atoms to preserve global fold
        # We'll add harmonic restraints on CA atoms to not stray too far (kcal/mol/Å^2 converted)
        from openmm import CustomExternalForce
        k_rest = 10.0 * unit.kilocalories_per_mole / (unit.angstrom**2)  # fairly strong restraint
        restraint = CustomExternalForce("0.5 * k * ((x - x0)^2 + (y - y0)^2 + (z - z0)^2)")
        restraint.addPerParticleParameter("x0")
        restraint.addPerParticleParameter("y0")
        restraint.addPerParticleParameter("z0")
        restraint.addGlobalParameter("k", k_rest)
        # add CA restraints
        # Build map atom index -> position
        for atom in topology.atoms():
            name = atom.name
            if name == "CA":
                idx = atom.index
                pos = positions[idx]
                restraint.addParticle(idx, [pos.x, pos.y, pos.z])
        system.addForce(restraint)

        # set up integrator and simulation
        integrator = LangevinIntegrator(300*unit.kelvin, 1.0/unit.picoseconds, 0.002*unit.picoseconds)
        platform = None
        if self.platform_name:
            platform = mm.Platform.getPlatformByName(self.platform_name)
        simulation = app.Simulation(topology, system, integrator, platform) if platform else app.Simulation(topology, system, integrator)
        simulation.context.setPositions(positions)

        # Minimize
        simulation.minimizeEnergy(maxIterations=niter_minimize)

        state = simulation.context.getState(getPositions=True, enforcePeriodicBox=False)
        minimized_positions = state.getPositions()

        # write out minimized pdb text
        out_pdb_path = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb").name
        with open(out_pdb_path, 'w') as f:
            app.PDBFile.writeFile(simulation.topology, minimized_positions, f)

        with open(out_pdb_path, 'r') as f:
            out_pdb_text = f.read()
        os.remove(out_pdb_path)
        return out_pdb_text

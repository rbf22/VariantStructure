from pdbfixer import PDBFixer
from openmm import app, unit
import openmm as mm
import tempfile, os
from Bio.PDB import PDBIO

def biopy_structure_to_pdb_text(structure):
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
    def __init__(self, ph=7.0, platform_name=None):
        self.ph = ph
        self.platform_name = platform_name

    def fix_and_minimize(self, biopy_struct, keep_water=True, niter_minimize=500, restrain_calpha=True):
        pdb_text = biopy_structure_to_pdb_text(biopy_struct)
        fixer = PDBFixer(pdbfile=None)
        fixer.readPDB(pdb_text)
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(self.ph)
        if not keep_water:
            fixer.removeHeterogens(keepWater=False)
        topology = fixer.topology
        positions = fixer.positions

        forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3p.xml')
        system = forcefield.createSystem(topology,
                                         nonbondedMethod=app.PME,
                                         nonbondedCutoff=1.0*unit.nanometer,
                                         constraints=app.HBonds)

        if restrain_calpha:
            from openmm import CustomExternalForce
            k_rest = 10.0 * unit.kilocalories_per_mole / (unit.angstrom**2)
            restraint = CustomExternalForce("0.5 * k * ((x - x0)^2 + (y - y0)^2 + (z - z0)^2)")
            restraint.addPerParticleParameter("x0")
            restraint.addPerParticleParameter("y0")
            restraint.addGlobalParameter("k", k_rest)
            for atom in topology.atoms():
                if atom.name == "CA":
                    idx = atom.index
                    pos = positions[idx]
                    restraint.addParticle(idx, [pos.x, pos.y, pos.z])
            system.addForce(restraint)

        integrator = app.LangevinIntegrator(300*unit.kelvin, 1.0/unit.picoseconds, 0.002*unit.picoseconds)
        platform = None
        if self.platform_name:
            platform = mm.Platform.getPlatformByName(self.platform_name)
        simulation = app.Simulation(topology, system, integrator, platform) if platform else app.Simulation(topology, system, integrator)
        simulation.context.setPositions(positions)

        simulation.minimizeEnergy(maxIterations=niter_minimize)
        state = simulation.context.getState(getPositions=True, enforcePeriodicBox=False)
        minimized_positions = state.getPositions()

        out_pdb_path = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb").name
        with open(out_pdb_path, 'w') as f:
            app.PDBFile.writeFile(simulation.topology, minimized_positions, f)

        with open(out_pdb_path, 'r') as f:
            out_pdb_text = f.read()
        os.remove(out_pdb_path)
        return out_pdb_text
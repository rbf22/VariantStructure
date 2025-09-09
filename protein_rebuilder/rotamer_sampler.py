from openmm import app, unit
import openmm as mm
import numpy as np
import random

class RotamerMC:
    def __init__(self, biopy_struct, dunbrack_lib, temperature=300):
        self.struct = biopy_struct
        self.dlib = dunbrack_lib
        self.temperature = temperature
        # Build OpenMM system & context once
        fixer = PDBFixer(pdbfile=None)
        fixer.readPDB(self._write_pdb_text(self.struct))
        fixer.findMissingAtoms(); fixer.addMissingAtoms(); fixer.addMissingHydrogens(7.0)
        ff = app.ForceField('amber14-all.xml', 'amber14/tip3p.xml')
        self.topo = fixer.topology
        self.positions = fixer.positions
        self.system = ff.createSystem(self.topo, nonbondedMethod=app.NoCutoff, constraints=None)
        integrator = mm.LangevinIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 0.002*unit.picoseconds)
        self.sim = app.Simulation(self.topo, self.system, integrator)
        self.sim.context.setPositions(self.positions)
        state = self.sim.context.getState(getEnergy=True)
        self.current_energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)

    def _write_pdb_text(self, struct):
        import tempfile
        from Bio.PDB import PDBIO
        fh = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb")
        fh.close()
        io = PDBIO()
        io.set_structure(struct)
        io.save(fh.name)
        txt = open(fh.name).read()
        import os
        os.remove(fh.name)
        return txt

    def sample(self, niter=500):
        kB = 0.0083144621  # kJ/mol/K
        model = self.struct[0]
        for it in range(niter):
            # pick a random residue
            # (chain, resid) selection omitted for brevity
            residue = ...
            # compute its backbone phi, psi (using e.g. Bio.PDB calc dihedrals)
            phi, psi = compute_phi_psi(residue)
            resname = residue.get_resname()
            rotamers = self.dlib.get_rotamers(resname, phi, psi)
            if not rotamers:
                continue
            choice = random.choices(rotamers, weights=[r['prob'] for r in rotamers])[0]
            # backup coords
            backup = {a.get_name(): a.get_coord().copy() for a in residue}
            # apply mean chi rotation
            set_chi_angles_for_residue(residue, choice['chi'])
            # evaluate energy
            state = self.sim.context.getState(getEnergy=True)
            E_new = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            dE = E_new - self.current_energy
            if dE <= 0 or random.random() < math.exp(-dE/(kB*self.temperature)):
                self.current_energy = E_new
            else:
                # revert
                for atom in residue:
                    atom.set_coord(backup[atom.get_name()])
        # at end, return final structure (self.struct)
        return self.struct
"""
A Monte Carlo sampler for protein side-chain rotamers.
"""

import math
import random
from io import StringIO
import tempfile
import os
from Bio.PDB import PDBIO, Polypeptide
from openmm import app, unit
import openmm as mm
from pdbfixer import PDBFixer


def set_chi_angles_for_residue(_, __):
    """
    Placeholder for setting chi angles.
    Needs to be implemented.
    """


class RotamerMC:  # pylint: disable=too-many-instance-attributes
    """
    A Monte Carlo sampler for protein side-chain rotamers.
    """

    def __init__(self, biopy_struct, dunbrack_lib, temperature=300):
        self.struct = biopy_struct
        self.dlib = dunbrack_lib
        self.temperature = temperature
        # Build OpenMM system & context once
        sio = StringIO(self._write_pdb_text(self.struct))
        fixer = PDBFixer(pdbfile=sio)
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.0)
        ff = app.ForceField("amber14-all.xml", "amber14/tip3p.xml")
        self.topo = fixer.topology
        self.positions = fixer.positions
        self.system = ff.createSystem(
            self.topo, nonbondedMethod=app.NoCutoff, constraints=None
        )
        integrator = mm.LangevinIntegrator(
            300 * unit.kelvin,
            1.0 / unit.picosecond,  # pylint: disable=no-member
            0.002 * unit.picoseconds,  # pylint: disable=no-member
        )
        self.sim = app.Simulation(self.topo, self.system, integrator)
        self.sim.context.setPositions(self.positions)
        state = self.sim.context.getState(getEnergy=True)
        self.current_energy = state.getPotentialEnergy().value_in_unit(
            unit.kilojoule_per_mole
        )

    def _write_pdb_text(self, struct):
        with tempfile.NamedTemporaryFile(
            delete=False, suffix=".pdb", mode="w", encoding="utf-8"
        ) as fh:
            io = PDBIO()
            io.set_structure(struct)
            io.save(fh.name)
            tmpname = fh.name
        with open(tmpname, "r", encoding="utf-8") as f:
            txt = f.read()
        os.remove(tmpname)
        return txt

    def sample(self, niter=500):  # pylint: disable=too-many-locals
        """
        Sample rotamers for the protein.
        """
        k_b = 0.0083144621  # kJ/mol/K
        chain = self.struct[0]["A"]
        poly = Polypeptide.Polypeptide(chain)
        phi_psi_list = poly.get_phi_psi_list()
        for _ in range(niter):
            # pick a random residue
            residues = [r for r in chain.get_residues() if r.get_resname() != "HOH"]
            residue_index, residue = random.choice(list(enumerate(residues)))

            phi, psi = phi_psi_list[residue_index]
            if phi is None or psi is None:
                continue

            phi_deg = math.degrees(phi)
            psi_deg = math.degrees(psi)

            resname = residue.get_resname()
            rotamers = self.dlib.get_rotamers(resname, phi_deg, psi_deg)
            if not rotamers:
                continue
            choice = random.choices(rotamers, weights=[r["prob"] for r in rotamers])[0]
            # backup coords
            backup = {a.get_name(): a.get_coord().copy() for a in residue}
            # apply mean chi rotation
            set_chi_angles_for_residue(residue, choice["chi"])
            # evaluate energy
            state = self.sim.context.getState(getEnergy=True)
            e_new = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            d_e = e_new - self.current_energy
            if d_e <= 0 or random.random() < math.exp(-d_e / (k_b * self.temperature)):
                self.current_energy = e_new
            else:
                # revert
                for atom in residue:
                    atom.set_coord(backup[atom.get_name()])
        # at end, return final structure (self.struct)
        return self.struct

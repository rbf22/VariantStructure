"""
A Monte Carlo sampler for protein side-chain rotamers.
"""

import math
import os
import random
import tempfile
from io import StringIO

import numpy as np
from Bio.PDB import PDBIO, Polypeptide
from Bio.PDB.vectors import rotaxis, calc_dihedral
import openmm as mm
from openmm import app, unit
from pdbfixer import PDBFixer

CHI_ATOMS = {
    # CHI1
    ("ARG", 1): ("N", "CA", "CB", "CG"),
    ("ASN", 1): ("N", "CA", "CB", "CG"),
    ("ASP", 1): ("N", "CA", "CB", "CG"),
    ("CYS", 1): ("N", "CA", "CB", "SG"),
    ("GLN", 1): ("N", "CA", "CB", "CG"),
    ("GLU", 1): ("N", "CA", "CB", "CG"),
    ("HIS", 1): ("N", "CA", "CB", "CG"),
    ("ILE", 1): ("N", "CA", "CB", "CG1"),
    ("LEU", 1): ("N", "CA", "CB", "CG"),
    ("LYS", 1): ("N", "CA", "CB", "CG"),
    ("MET", 1): ("N", "CA", "CB", "CG"),
    ("PHE", 1): ("N", "CA", "CB", "CG"),
    ("PRO", 1): ("N", "CA", "CB", "CG"),
    ("SER", 1): ("N", "CA", "CB", "OG"),
    ("THR", 1): ("N", "CA", "CB", "OG1"),
    ("TRP", 1): ("N", "CA", "CB", "CG"),
    ("TYR", 1): ("N", "CA", "CB", "CG"),
    ("VAL", 1): ("N", "CA", "CB", "CG1"),
    # CHI2
    ("ARG", 2): ("CA", "CB", "CG", "CD"),
    ("ASN", 2): ("CA", "CB", "CG", "OD1"),
    ("ASP", 2): ("CA", "CB", "CG", "OD1"),
    ("GLN", 2): ("CA", "CB", "CG", "CD"),
    ("GLU", 2): ("CA", "CB", "CG", "CD"),
    ("HIS", 2): ("CA", "CB", "CG", "ND1"),
    ("ILE", 2): ("CA", "CB", "CG1", "CD1"),
    ("LEU", 2): ("CA", "CB", "CG", "CD1"),
    ("LYS", 2): ("CA", "CB", "CG", "CD"),
    ("MET", 2): ("CA", "CB", "CG", "SD"),
    ("PHE", 2): ("CA", "CB", "CG", "CD1"),
    ("PRO", 2): ("CA", "CB", "CG", "CD"),
    ("TRP", 2): ("CA", "CB", "CG", "CD1"),
    ("TYR", 2): ("CA", "CB", "CG", "CD1"),
    # CHI3
    ("ARG", 3): ("CB", "CG", "CD", "NE"),
    ("GLN", 3): ("CB", "CG", "CD", "OE1"),
    ("GLU", 3): ("CB", "CG", "CD", "OE1"),
    ("LYS", 3): ("CB", "CG", "CD", "CE"),
    ("MET", 3): ("CB", "CG", "SD", "CE"),
    # CHI4
    ("ARG", 4): ("CG", "CD", "NE", "CZ"),
    ("LYS", 4): ("CG", "CD", "CE", "NZ"),
    # CHI5
    ("ARG", 5): ("CD", "NE", "CZ", "NH1"),
}

SIDECHAIN_TREE = {
    "CYS": {"CB": ["SG"]},
    "ASP": {"CB": ["CG"], "CG": ["OD1", "OD2"]},
    "SER": {"CB": ["OG"]},
    "GLN": {"CB": ["CG"], "CG": ["CD"], "CD": ["OE1", "NE2"]},
    "LYS": {"CB": ["CG"], "CG": ["CD"], "CD": ["CE"], "CE": ["NZ"]},
    "ILE": {"CB": ["CG1", "CG2"], "CG1": ["CD1"]},
    "PRO": {"CB": ["CG"], "CG": ["CD"]},
    "THR": {"CB": ["OG1", "CG2"]},
    "PHE": {
        "CB": ["CG"],
        "CG": ["CD1", "CD2"],
        "CD1": ["CE1"],
        "CD2": ["CE2"],
        "CE1": ["CZ"],
        "CE2": ["CZ"],
    },
    "ASN": {"CB": ["CG"], "CG": ["OD1", "ND2"]},
    "GLY": {},
    "HIS": {
        "CB": ["CG"],
        "CG": ["ND1", "CD2"],
        "ND1": ["CE1"],
        "CD2": ["NE2"],
        "CE1": ["NE2"],
    },
    "LEU": {"CB": ["CG"], "CG": ["CD1", "CD2"]},
    "ARG": {
        "CB": ["CG"],
        "CG": ["CD"],
        "CD": ["NE"],
        "NE": ["CZ"],
        "CZ": ["NH1", "NH2"],
    },
    "TRP": {
        "CB": ["CG"],
        "CG": ["CD1", "CD2"],
        "CD1": ["NE1"],
        "CD2": ["CE2"],
        "NE1": ["CE2"],
        "CE2": ["CZ2", "CH2", "CZ3"],
        "CZ2": ["CH2"],
        "CH2": ["CZ3"],
    },
    "VAL": {"CB": ["CG1", "CG2"]},
    "GLU": {"CB": ["CG"], "CG": ["CD"], "CD": ["OE1", "OE2"]},
    "MET": {"CB": ["CG"], "CG": ["SD"], "SD": ["CE"]},
    "TYR": {
        "CB": ["CG"],
        "CG": ["CD1", "CD2"],
        "CD1": ["CE1"],
        "CD2": ["CE2"],
        "CE1": ["CZ"],
        "CE2": ["CZ"],
        "CZ": ["OH"],
    },
}


def get_downstream_atoms(residue, start_atom):
    """Get all atoms downstream of a starting atom in the sidechain."""
    resname = residue.get_resname()
    if resname not in SIDECHAIN_TREE:
        return []

    downstream = set()
    q = [start_atom]
    visited = {start_atom}

    while q:
        curr = q.pop(0)
        if curr in SIDECHAIN_TREE[resname]:
            for child in SIDECHAIN_TREE[resname][curr]:
                if child not in visited:
                    downstream.add(child)
                    q.append(child)
                    visited.add(child)

    return [residue[atom_name] for atom_name in downstream if atom_name in residue]


def set_chi_angles_for_residue(residue, chi_angles):
    """
    Set the chi angles for a residue.
    """
    resname = residue.get_resname()
    for i, angle_deg in enumerate(chi_angles):
        chi_index = i + 1
        key = (resname, chi_index)
        if key not in CHI_ATOMS:
            continue

        a, b, c, d = CHI_ATOMS[key]

        try:
            v_a = residue[a].get_vector()
            v_b = residue[b].get_vector()
            v_c = residue[c].get_vector()
            v_d = residue[d].get_vector()
        except KeyError:
            # Atom not found, cannot set angle
            continue

        current_angle_rad = calc_dihedral(v_a, v_b, v_c, v_d)
        current_angle_deg = np.degrees(current_angle_rad)

        # The angle from the rotamer library is what we want to set it to.
        target_angle_deg = angle_deg

        rotation_angle_deg = target_angle_deg - current_angle_deg
        rotation_angle_rad = np.radians(rotation_angle_deg)

        rotation_axis = (v_c - v_b).normalized()

        # The atoms to rotate are d and everything downstream of it.
        atoms_to_rotate = get_downstream_atoms(residue, d)
        if d in residue:
            atoms_to_rotate.insert(0, residue[d])

        rot_matrix = rotaxis(rotation_angle_rad, rotation_axis)

        for atom in atoms_to_rotate:
            atom_vec = atom.get_vector()
            new_pos = (atom_vec - v_b).left_multiply(rot_matrix) + v_b
            atom.set_coord(new_pos.get_array())


class RotamerMC:  # pylint: disable=too-many-instance-attributes,too-few-public-methods
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

    def sample(self, niter=500):
        """
        Sample rotamers for the protein.
        """
        # pylint: disable=too-many-locals
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

            # Update positions in OpenMM context
            self._update_omm_positions(residue)

            # evaluate energy
            state = self.sim.context.getState(getEnergy=True)
            e_new = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            d_e = e_new - self.current_energy

            if d_e <= 0 or random.random() < math.exp(-d_e / (k_b * self.temperature)):
                self.current_energy = e_new
                # Keep the new positions in self.positions
                self.positions = self.sim.context.getState(
                    getPositions=True
                ).getPositions(asNumpy=True)
            else:
                # revert
                for atom in residue:
                    atom.set_coord(backup[atom.get_name()])
                self._update_omm_positions(residue)

        # at end, return final structure (self.struct)
        final_pos = self.sim.context.getState(getPositions=True).getPositions()

        # To be safe, update the biopython structure with the final positions from OMM
        for atom in self.struct.get_atoms():
            res_id = atom.get_parent().get_id()
            atom_name = atom.get_name()

            # Find corresponding atom in OMM topology
            for res_topo in self.topo.residues():
                if (
                    res_topo.name == atom.get_parent().get_resname()
                    and res_topo.id == str(res_id[1])
                ):
                    for atom_topo in res_topo.atoms():
                        if atom_topo.name == atom_name:
                            atom.set_coord(
                                final_pos[atom_topo.index].value_in_unit(unit.angstrom)
                            )
                            break
                    else:
                        continue
                    break
        return self.struct

    def _update_omm_positions(self, residue):
        """
        Update the OpenMM positions from a Biopython residue.
        """
        res_id = residue.get_id()

        # Create a copy of the current positions
        new_positions = self.positions.value_in_unit(
            unit.nanometer  # pylint: disable=no-member
        ).copy()

        # Find the residue in the topology
        for res_topo in self.topo.residues():
            if res_topo.name == residue.get_resname() and res_topo.id == str(res_id[1]):
                for atom_topo in res_topo.atoms():
                    if atom_topo.name in residue:
                        new_coord = residue[atom_topo.name].get_coord()
                        new_positions[atom_topo.index] = (
                            new_coord / 10.0
                        )  # Angstrom to nm

        self.sim.context.setPositions(new_positions)

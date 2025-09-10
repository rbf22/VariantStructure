"""Helper functions for reading and writing PDB files."""

from Bio.PDB import PDBParser, PDBIO, Select


def read_pdb_file(pdb_path_or_handle):
    """Read a PDB file and return a Bio.PDB.Structure object."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("ref", pdb_path_or_handle)
    return structure


class _KeepAll(Select):
    """A PDBIO Select class that keeps all residues."""

    def accept_residue(self, residue):
        """Accept all residues."""
        return True


def write_pdb_file(structure, out_path):
    """Write a Bio.PDB.Structure object to a PDB file."""
    io_writer = PDBIO()
    io_writer.set_structure(structure)
    io_writer.save(out_path, select=_KeepAll())

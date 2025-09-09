from Bio.PDB import PDBParser, PDBIO, Select

def read_pdb_file(pdb_path_or_handle):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('ref', pdb_path_or_handle)
    return structure

class _KeepAll(Select):
    def accept_residue(self, residue): return True

def write_pdb_file(structure, out_path):
    io_writer = PDBIO()
    io_writer.set_structure(structure)
    io_writer.save(out_path, select=_KeepAll())
import os
from protein_rebuilder.io_helpers import read_pdb_file, write_pdb_file

def test_write_pdb_file():
    """
    Test that a PDB file can be written.
    """
    test_pdb = os.path.join(os.path.dirname(__file__), "data", "4EFL_snippet.pdb")
    struct = read_pdb_file(test_pdb)
    output_pdb = os.path.join(os.path.dirname(__file__), "data", "output.pdb")
    write_pdb_file(struct, output_pdb)
    assert os.path.exists(output_pdb)
    # Read the file back and check that it is the same
    struct2 = read_pdb_file(output_pdb)
    assert len(list(struct.get_atoms())) == len(list(struct2.get_atoms()))
    os.remove(output_pdb)

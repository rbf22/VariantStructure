import os
from protein_rebuilder.io_helpers import read_pdb_file
from protein_rebuilder.rotamer_sampler import RotamerMC
from protein_rebuilder.rotamer_loader import DunbrackLibrary


def test_rotamer_mc():
    """
    Test the RotamerMC class.
    """
    test_pdb = os.path.join(os.path.dirname(__file__), "data", "4EFL_snippet.pdb")
    struct = read_pdb_file(test_pdb)
    dunbrack_lib = DunbrackLibrary()
    rotamer_mc = RotamerMC(struct, dunbrack_lib)
    assert rotamer_mc is not None
    final_struct = rotamer_mc.sample(niter=10)
    assert final_struct is not None

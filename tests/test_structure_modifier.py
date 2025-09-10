import os
from protein_rebuilder.io_helpers import read_pdb_file
from protein_rebuilder.structure_modifier import get_chain_sequence


def test_get_chain_sequence():
    """
    Test that the sequence of a chain can be retrieved.
    """
    test_pdb = os.path.join(os.path.dirname(__file__), "data", "4EFL_snippet.pdb")
    struct = read_pdb_file(test_pdb)
    chain = struct[0]["A"]
    seq, _ = get_chain_sequence(chain)
    assert seq == "MTEYK"

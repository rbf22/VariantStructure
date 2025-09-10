import numpy as np
from protein_rebuilder.ccd_loop_builder import build_peptide_fragment, length


def test_build_peptide_fragment():
    """
    Test that a peptide fragment can be built.
    """
    n_residues = 3
    fragment = build_peptide_fragment(n_residues)
    assert len(fragment) == n_residues
    # Check that the distance between CA atoms is roughly correct
    for i in range(n_residues - 1):
        ca1 = fragment[i]["CA"]
        ca2 = fragment[i + 1]["CA"]
        dist = length(ca2 - ca1)
        assert np.isclose(dist, 3.8, atol=0.1)

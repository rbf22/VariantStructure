from protein_rebuilder.rotamer_lib import get_rotamers_for_resname


def test_get_rotamers_for_resname():
    """
    Test that rotamers can be retrieved for a residue name.
    """
    rotamers = get_rotamers_for_resname("ALA")
    assert isinstance(rotamers, list)

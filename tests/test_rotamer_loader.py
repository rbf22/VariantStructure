from protein_rebuilder.rotamer_loader import DunbrackLibrary

def test_dunbrack_library():
    """
    Test that the Dunbrack library can be loaded.
    """
    lib = DunbrackLibrary()
    assert lib is not None
    # Check that we can get rotamers for a known residue
    rotamers = lib.get_rotamers("VAL", -60, 140)
    assert isinstance(rotamers, list)
    assert len(rotamers) > 0

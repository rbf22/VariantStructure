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


def test_dunbrack_library_parsing():
    """
    Test that the Dunbrack library is parsed correctly.
    """
    lib = DunbrackLibrary()
    # Check the first rotamer for ARG at phi=-180, psi=-180
    rotamers = lib.get_rotamers("ARG", -180, -180)
    assert len(rotamers) > 0
    assert rotamers[0]["prob"] == 0.24973


def test_all_amino_acids_have_rotamers():
    """
    Test that all standard amino acids have rotamers in the library.
    """
    lib = DunbrackLibrary()
    amino_acids = [
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
    ]
    for aa in amino_acids:
        rotamers = lib.get_rotamers(aa, -60, 140)
        assert isinstance(rotamers, list)
        assert len(rotamers) > 0

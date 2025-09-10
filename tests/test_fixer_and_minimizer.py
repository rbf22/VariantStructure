import os
from protein_rebuilder.fixer_and_minimizer import biopy_structure_to_pdb_text
from protein_rebuilder.io_helpers import read_pdb_file
from protein_rebuilder.fixer_and_minimizer import FixerMinimizer
from protein_rebuilder.structure_modifier import StructureModifier, get_chain_sequence
from protein_rebuilder.sequence_aligner import SequenceAligner


def test_biopy_structure_to_pdb_text():
    """
    Test that a Biopython structure can be converted to a PDB file string.
    """
    test_pdb = os.path.join(os.path.dirname(__file__), "data", "4EFL_snippet.pdb")
    struct = read_pdb_file(test_pdb)
    pdb_text = biopy_structure_to_pdb_text(struct)
    assert "ATOM" in pdb_text
    assert "TER" in pdb_text
    assert "END" in pdb_text


def test_fixer_minimizer():
    """
    Test the FixerMinimizer class.
    """
    test_pdb = os.path.join(os.path.dirname(__file__), "data", "4EFL_snippet.pdb")
    struct = read_pdb_file(test_pdb)
    chain = struct[0]["A"]
    ref_seq, _ = get_chain_sequence(chain)
    # create new seq with a small insertion
    new_seq = ref_seq[:2] + "GG" + ref_seq[2:]
    a_ref, a_new = SequenceAligner(ref_seq, new_seq).align()
    mod = StructureModifier(struct, chain_id="A")
    new_struct = mod.build_modified_structure(a_ref, a_new)
    fm = FixerMinimizer()
    pdb_text = fm.fix_and_minimize(new_struct)
    assert "ATOM" in pdb_text

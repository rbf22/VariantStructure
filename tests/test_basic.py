import os
from protein_rebuilder.io_helpers import read_pdb_file
from protein_rebuilder.sequence_aligner import SequenceAligner
from protein_rebuilder.structure_modifier import StructureModifier, get_chain_sequence


def test_read_pdb():
    """
    Test that a PDB file can be read.
    """
    test_pdb = os.path.join(os.path.dirname(__file__), "data", "4EFL_snippet.pdb")
    struct = read_pdb_file(test_pdb)
    assert struct is not None


def test_alignment_simple():
    """
    Test a simple sequence alignment.
    """
    ref = "MKT"
    new = "MKT"
    a_ref, a_new = SequenceAligner(ref, new).align()
    assert a_ref.replace("-", "") == ref


def test_structure_modifier_smoke():
    """
    Test that the structure modifier runs without errors.
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
    assert new_struct is not None


def test_cli():
    """
    Test the command-line interface.
    """
    test_pdb = os.path.join(os.path.dirname(__file__), "data", "4EFL_snippet.pdb")
    new_seq_file = os.path.join(os.path.dirname(__file__), "data", "new_seq.fasta")
    with open(new_seq_file, "w") as f:
        f.write(">new_seq\n")
        f.write("MTEYKL")

    output_pdb = os.path.join(os.path.dirname(__file__), "data", "output.pdb")

    from protein_rebuilder.cli import main

    main(
        [
            "--pdb",
            test_pdb,
            "--new-seq",
            new_seq_file,
            "--out",
            output_pdb,
        ]
    )
    assert os.path.exists(output_pdb)

    # Clean up the output file
    os.remove(new_seq_file)
    os.remove(output_pdb)

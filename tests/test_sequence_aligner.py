from protein_rebuilder.sequence_aligner import SequenceAligner


def test_sequence_aligner():
    """
    Test the SequenceAligner class.
    """
    ref_seq = "MTEYK"
    new_seq = "MTEYKL"
    aligner = SequenceAligner(ref_seq, new_seq)
    a_ref, a_new = aligner.align()
    assert a_ref == "MTEYK-"
    assert a_new == "MTEYKL"

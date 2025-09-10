"""
A simple sequence aligner using Biopython's pairwise2 module.
"""

from Bio import pairwise2


class SequenceAligner:  # pylint: disable=too-few-public-methods
    """
    A class to align two sequences.
    """

    def __init__(self, ref_seq, new_seq):
        self.ref_seq = ref_seq
        self.new_seq = new_seq

    def align(self):
        """
        Align the two sequences.
        """
        aln = pairwise2.align.globalxx(
            self.ref_seq, self.new_seq, one_alignment_only=True
        )
        if len(aln) == 0:
            raise RuntimeError("No alignment returned")
        a_ref, a_new, _, _, _ = aln[0]
        return a_ref, a_new

# protein_rebuilder/sequence_aligner.py
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.SeqUtils import seq1

class SequenceAligner:
    """
    Aligns reference sequence (from PDB chain) to new sequence and
    produces aligned pairs for mapping.
    """

    def __init__(self, ref_seq, new_seq):
        """
        ref_seq: string, single-letter AA sequence from reference structure
        new_seq: string, single-letter AA sequence to model
        """
        self.ref_seq = ref_seq
        self.new_seq = new_seq

    def align(self):
        # Use a simple global alignment (identity scoring); adapt scoring if desired.
        aln = pairwise2.align.globalxx(self.ref_seq, self.new_seq, one_alignment_only=True)
        if len(aln) == 0:
            raise RuntimeError("No alignment returned")
        a_ref, a_new, score, begin, end = aln[0]
        return a_ref, a_new

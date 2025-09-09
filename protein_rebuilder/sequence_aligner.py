from Bio import pairwise2

class SequenceAligner:
    def __init__(self, ref_seq, new_seq):
        self.ref_seq = ref_seq
        self.new_seq = new_seq

    def align(self):
        aln = pairwise2.align.globalxx(self.ref_seq, self.new_seq, one_alignment_only=True)
        if len(aln) == 0:
            raise RuntimeError("No alignment returned")
        a_ref, a_new, score, begin, end = aln[0]
        return a_ref, a_new
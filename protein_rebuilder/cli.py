"""
Command-line interface for protein_rebuilder.
"""

import argparse
from Bio.PDB import PDBParser
from .io_helpers import read_pdb_file, write_pdb_file
from .sequence_aligner import SequenceAligner
from .structure_modifier import StructureModifier, get_chain_sequence
from .fixer_and_minimizer import FixerMinimizer
from .rotamer_sampler import RotamerMC
from .rotamer_loader import DunbrackLibrary


def parse_args(argv=None):
    """
    Parse command-line arguments.
    """
    p = argparse.ArgumentParser(prog="protein_rebuilder")
    p.add_argument("rebuild", nargs="?", help="subcommand", default="rebuild")
    p.add_argument("--pdb", required=True, help="reference pdb file")
    p.add_argument("--chain", default="A", help="chain to operate on")
    p.add_argument(
        "--new-seq",
        required=True,
        help="path to FASTA or plain sequence file (single sequence)",
    )
    p.add_argument("--out", default="rebuilt.pdb", help="output pdb file")
    p.add_argument(
        "--min-iterations", type=int, default=500, help="OpenMM minimization iterations"
    )
    p.add_argument(
        "--mc-iterations", type=int, default=200, help="rotamer Monte Carlo iterations"
    )
    return p.parse_args(argv)


def read_seq_from_file(path):
    """
    Read a sequence from a file.
    """
    with open(path, "r", encoding="utf-8") as f:
        data = f.read().strip()
    if "\n" in data:
        # assume FASTA: take first non-header line(s)
        lines = [line.strip() for line in data.splitlines() if not line.startswith(">")]
        return "".join(lines)
    return data


def _run_alignment_and_modification(struct, new_seq, chain_id):
    """Align sequences and build the initial modified structure."""
    chain = struct[0][chain_id]
    ref_seq, _ = get_chain_sequence(chain)
    print("Ref seq length:", len(ref_seq))
    aligner = SequenceAligner(ref_seq, new_seq)
    a_ref, a_new = aligner.align()
    modifier = StructureModifier(struct, chain_id=chain_id)
    return modifier.build_modified_structure(a_ref, a_new)


def _run_minimization(modified_struct, min_iterations):
    """Fix and minimize a structure, returning PDB text."""
    fm = FixerMinimizer()
    print("Running initial fix & minimize...")
    pdb_text = fm.fix_and_minimize(
        modified_struct, keep_water=False, niter_minimize=min_iterations
    )
    with open("intermediate_fixed.pdb", "w", encoding="utf-8") as fh:
        fh.write(pdb_text)
    return pdb_text


def _run_rotamer_sampling(mc_iterations):
    """Run rotamer sampling and return the final structure."""
    print("Running rotamer Monte Carlo...")
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("fixed", "intermediate_fixed.pdb")
    dunbrack_lib = DunbrackLibrary()
    rotamer_mc = RotamerMC(struct, dunbrack_lib)
    return rotamer_mc.sample(niter=mc_iterations)


def main(argv=None):
    """
    Main function for the command-line interface.
    """
    args = parse_args(argv)

    print("Reading PDB:", args.pdb)
    struct = read_pdb_file(args.pdb)
    new_seq = read_seq_from_file(args.new_seq)

    modified = _run_alignment_and_modification(struct, new_seq, args.chain)

    _run_minimization(modified, args.min_iterations)

    final_struct = _run_rotamer_sampling(args.mc_iterations)

    write_pdb_file(final_struct, args.out)
    print("Wrote:", args.out)

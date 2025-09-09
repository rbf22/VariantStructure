import argparse
from .io_helpers import read_pdb_file, write_pdb_file
from .sequence_aligner import SequenceAligner
from .structure_modifier import StructureModifier
from .fixer_and_minimizer import FixerMinimizer
from .rotamer_sampler import rotamer_monte_carlo
import sys

def parse_args(argv=None):
    p = argparse.ArgumentParser(prog="protein_rebuilder")
    p.add_argument("rebuild", nargs='?', help="subcommand", default="rebuild")
    p.add_argument("--pdb", required=True, help="reference pdb file")
    p.add_argument("--chain", default="A", help="chain to operate on")
    p.add_argument("--new-seq", required=True, help="path to FASTA or plain sequence file (single sequence)")
    p.add_argument("--out", default="rebuilt.pdb", help="output pdb file")
    p.add_argument("--min-iterations", type=int, default=500, help="OpenMM minimization iterations")
    p.add_argument("--mc-iterations", type=int, default=200, help="rotamer Monte Carlo iterations")
    return p.parse_args(argv)

def read_seq_from_file(path):
    with open(path) as f:
        data = f.read().strip()
    if '\n' in data:
        # assume FASTA: take first non-header line(s)
        lines = [l.strip() for l in data.splitlines() if not l.startswith('>')]
        return ''.join(lines)
    return data

def main(argv=None):
    args = parse_args(argv)
    pdb_file = args.pdb
    new_seq = read_seq_from_file(args.new_seq)
    print("Reading PDB:", pdb_file)
    struct = read_pdb_file(pdb_file)
    chain = struct[0][args.chain]
    from .structure_modifier import get_chain_sequence
    ref_seq, _ = get_chain_sequence(chain)
    print("Ref seq length:", len(ref_seq))
    from .sequence_aligner import SequenceAligner
    aligner = SequenceAligner(ref_seq, new_seq)
    a_ref, a_new = aligner.align()
    modifier = StructureModifier(struct, chain_id=args.chain)
    modified = modifier.build_modified_structure(a_ref, a_new)
    # run PDBFixer+OpenMM minimization to add missing atoms
    fm = FixerMinimizer()
    print("Running initial fix & minimize...")
    pdb_text = fm.fix_and_minimize(modified, keep_water=False, niter_minimize=args.min_iterations)
    # Write intermediate
    with open("intermediate_fixed.pdb","w") as fh:
        fh.write(pdb_text)
    # Rotamer Monte Carlo
    print("Running rotamer Monte Carlo...")
    # convert pdb_text back to biopy structure for rotamer sampler convenience
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    struct2 = parser.get_structure("fixed","intermediate_fixed.pdb")
    final_pdb_text = rotamer_monte_carlo(struct2, niter=args.mc_iterations)
    with open(args.out, "w") as fh:
        fh.write(final_pdb_text)
    print("Wrote:", args.out)
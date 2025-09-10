# protein_rebuilder/example_hras.py
"""
Example usage:
- Downloads a PDB of H-Ras (4EFL)
- Reads chain A sequence
- Aligns a provided new sequence (example: small mutation + insertion)
- Produces a rebuilt PDB with sidechains and minimized structure
"""

import urllib.request
from protein_rebuilder.sequence_aligner import SequenceAligner
from protein_rebuilder.structure_modifier import StructureModifier, get_chain_sequence
from protein_rebuilder.fixer_and_minimizer import FixerMinimizer
from protein_rebuilder.io_helpers import write_pdb_file, read_pdb_file


def download_pdb(pdb_id, outpath=None):
    pdb_id = pdb_id.lower()
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    if outpath is None:
        outpath = f"{pdb_id}.pdb"
    urllib.request.urlretrieve(url, outpath)
    return outpath


def run_example():
    pdb_id = "4EFL"  # H-Ras GppNHp-bound (example)
    pdb_file = download_pdb(pdb_id)
    print("Downloaded:", pdb_file)

    # Load structure
    struct = read_pdb_file(pdb_file)
    chain_id = "A"
    model = struct[0]
    chain = model[chain_id]

    # Get reference sequence
    ref_seq, _ = get_chain_sequence(chain)
    print("Reference sequence length:", len(ref_seq))

    # Example: create a new sequence with:
    # - mutation at position 12: G->V (common oncogenic)
    # - insertion after position 61: insert "GGG"
    # We'll produce a new sequence string for demonstration.
    # First create a copy of ref_seq and modify:
    # Note: positions are 1-based in biology; our string indices are 0-based.
    new_chars = list(ref_seq)
    # mutate G12 -> V
    if len(new_chars) >= 12:
        new_chars[11] = "V"
    # insert GGG after residue 61 if possible
    if len(new_chars) >= 61:
        new_chars = new_chars[:61] + list("GGG") + new_chars[61:]
    new_seq = "".join(new_chars)
    print("New sequence length:", len(new_seq))

    # Align
    aligner = SequenceAligner(ref_seq, new_seq)
    a_ref, a_new = aligner.align()
    print("Alignment length:", len(a_ref))

    # Build modified structure (placeholders for insertions)
    modifier = StructureModifier(struct, chain_id=chain_id)
    modified_struct = modifier.build_modified_structure(a_ref, a_new)

    # Save intermediate PDB (unfixed)
    write_pdb_file(modified_struct, "hras_modified_unfixed.pdb")
    print("Wrote intermediate unfixed PDB: hras_modified_unfixed.pdb")

    # Fix and minimize
    fm = FixerMinimizer(ph=7.0)
    minimized_pdb_text = fm.fix_and_minimize(
        modified_struct, keep_water=False, niter_minimize=500
    )
    # Save final PDB text
    with open("hras_rebuilt_minimized.pdb", "w") as f:
        f.write(minimized_pdb_text)
    print("Wrote rebuilt, minimized PDB: hras_rebuilt_minimized.pdb")
    print("Done.")


if __name__ == "__main__":
    run_example()

# protein_rebuilder/example_kinase.py
"""
Example usage with a kinase (Wee1) from PDB 9D0R.
- Downloads a PDB of Wee1 (9D0R)
- Reads chain A sequence
- Demonstrates various modifications and rebuilding.
"""

import urllib.request
from Bio.PDB import MMCIFParser, PDBParser
from protein_rebuilder.sequence_aligner import SequenceAligner
from protein_rebuilder.structure_modifier import StructureModifier, get_chain_sequence
from protein_rebuilder.io_helpers import write_pdb_file


def download_structure(pdb_id, outpath=None):
    """Downloads a structure file from the RCSB."""
    pdb_id = pdb_id.upper()

    # Try CIF first
    url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    if outpath is None:
        outpath = f"{pdb_id}.cif"

    try:
        urllib.request.urlretrieve(url, outpath)
        return outpath, "cif"
    except urllib.error.HTTPError as e:
        if e.code == 404:
            print(f"CIF format for {pdb_id} not found, trying PDB format.")
            url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
            pdb_outpath = f"{pdb_id}.pdb"
            urllib.request.urlretrieve(url, pdb_outpath)
            return pdb_outpath, "pdb"
        else:
            raise e


def read_structure_file(filepath, file_format):
    """Read a structure file (PDB or CIF) and return a Bio.PDB.Structure object."""
    if file_format == "cif":
        parser = MMCIFParser()
    else:
        parser = PDBParser(QUIET=True)

    structure = parser.get_structure("ref", filepath)
    return structure


def run_modification(struct, chain_id, ref_seq, new_seq, output_filename):
    """
    Runs the protein rebuilding process for a given modification.
    """
    print(f"\n--- Running modification for {output_filename} ---")
    print(f"Reference sequence length: {len(ref_seq)}")
    print(f"New sequence length: {len(new_seq)}")

    # Align sequences
    aligner = SequenceAligner(ref_seq, new_seq)
    a_ref, a_new = aligner.align()

    # Build modified structure
    modifier = StructureModifier(struct.copy(), chain_id=chain_id)
    modified_struct = modifier.build_modified_structure(a_ref, a_new)

    # Save the rebuilt PDB file
    write_pdb_file(modified_struct, output_filename)
    print(f"Wrote rebuilt PDB: {output_filename}")


def main():
    """
    Main function to run all kinase example modifications.
    """
    pdb_id = "9D0R"  # This is a Wee1 kinase
    struct_file, file_format = download_structure(pdb_id)
    print("Downloaded:", struct_file)

    # Load the base structure once
    struct = read_structure_file(struct_file, file_format)
    chain_id = "A"
    chain = struct[0][chain_id]

    # Get reference sequence and a map from residue number to sequence index
    ref_seq, res_list = get_chain_sequence(chain)
    resseq_to_idx = {res.id[1]: i for i, res in enumerate(res_list)}

    # --- Example 1: Missense alteration (e.g., T450A) ---
    res_to_mutate = 450
    idx_to_mutate = resseq_to_idx.get(res_to_mutate)
    if (
        idx_to_mutate is not None
        and len(ref_seq) > idx_to_mutate
        and ref_seq[idx_to_mutate] == "T"
    ):
        new_chars = list(ref_seq)
        new_chars[idx_to_mutate] = "A"
        new_seq_missense = "".join(new_chars)
        run_modification(
            struct, chain_id, ref_seq, new_seq_missense, "9d0r_missense_T450A.pdb"
        )

    # --- Example 2: Insertion (insert 'AGG' after residue 400) ---
    res_to_insert_after = 400
    idx_to_insert_after = resseq_to_idx.get(res_to_insert_after)
    if idx_to_insert_after is not None:
        new_chars = list(ref_seq)
        new_chars = (
            new_chars[: idx_to_insert_after + 1]
            + list("AGG")
            + new_chars[idx_to_insert_after + 1 :]
        )
        new_seq_insertion = "".join(new_chars)
        run_modification(
            struct, chain_id, ref_seq, new_seq_insertion, "9d0r_insertion_400AGG.pdb"
        )

    # --- Example 3: Deletion (delete residues 350-355) ---
    res_del_start = 350
    res_del_end = 355
    idx_del_start = resseq_to_idx.get(res_del_start)
    idx_del_end = resseq_to_idx.get(res_del_end)
    if idx_del_start is not None and idx_del_end is not None:
        new_chars = list(ref_seq)
        del new_chars[idx_del_start : idx_del_end + 1]
        new_seq_deletion = "".join(new_chars)
        run_modification(
            struct, chain_id, ref_seq, new_seq_deletion, "9d0r_deletion_350-355.pdb"
        )

    # --- Example 4: N-terminal truncation (remove first 10 residues) ---
    new_seq_n_trunc = ref_seq[10:]
    run_modification(struct, chain_id, ref_seq, new_seq_n_trunc, "9d0r_n_trunc_10.pdb")

    # --- Example 5: C-terminal truncation (remove last 10 residues) ---
    new_seq_c_trunc = ref_seq[:-10]
    run_modification(struct, chain_id, ref_seq, new_seq_c_trunc, "9d0r_c_trunc_10.pdb")

    # --- Example 6: Combination of alterations ---
    # T450A mutation, deletion of 350-355, and C-terminal truncation of 5 residues
    if (
        idx_to_mutate is not None
        and len(ref_seq) > idx_to_mutate
        and ref_seq[idx_to_mutate] == "T"
        and idx_del_start is not None
        and idx_del_end is not None
    ):
        new_chars = list(ref_seq)
        del new_chars[idx_del_start : idx_del_end + 1]

        idx_mut_shifted = resseq_to_idx.get(res_to_mutate)
        if idx_mut_shifted is not None and idx_mut_shifted > idx_del_end:
            idx_mut_shifted = idx_mut_shifted - (idx_del_end - idx_del_start + 1)

        if (
            idx_mut_shifted is not None
            and len(new_chars) > idx_mut_shifted
            and new_chars[idx_mut_shifted] == "T"
        ):
            new_chars[idx_mut_shifted] = "A"

        new_seq_combo = "".join(new_chars)
        new_seq_combo = new_seq_combo[:-5]
        run_modification(
            struct, chain_id, ref_seq, new_seq_combo, "9d0r_combination.pdb"
        )


if __name__ == "__main__":
    main()

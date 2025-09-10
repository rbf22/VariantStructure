"""
Modify a protein structure based on a new sequence.
"""

from Bio.Data.PDBData import protein_letters_3to1, protein_letters_1to3
from Bio.PDB import Model, Chain, Residue, Atom, Structure
import numpy as np
from .ccd_loop_builder import build_peptide_fragment, residues_to_biopy, ccd_close_loop


def get_chain_sequence(chain):
    """
    Get the amino acid sequence of a chain.
    """
    seq = []
    res_list = []
    for res in chain:
        try:
            aa = protein_letters_3to1[res.get_resname()]
        except KeyError:
            aa = "X"
        seq.append(aa)
        res_list.append(res)
    return "".join(seq), res_list


class StructureModifier:  # pylint: disable=too-few-public-methods
    """
    Modifies a protein structure based on a new sequence.
    """

    def __init__(self, structure, chain_id="A"):
        self.structure = structure
        self.chain_id = chain_id
        self.ca_positions = []

    def _clone_and_mutate_residue(self, new_chain, orig_res, new_res_id, new_aa):
        """Clone a residue and optionally mutate its residue name."""
        rid = (" ", new_res_id, " ")
        new_res = Residue.Residue(rid, orig_res.get_resname(), orig_res.get_segid())
        for atom in orig_res:
            at = Atom.Atom(
                atom.get_name(),
                atom.get_vector().get_array(),
                atom.get_bfactor(),
                atom.get_occupancy(),
                atom.get_altloc(),
                atom.get_fullname(),
                atom.get_serial_number(),
                element=atom.element,
            )
            new_res.add(at)
        # rename if needed
        if orig_res.get_resname() != protein_letters_1to3.get(new_aa):
            try:
                new_res.resname = protein_letters_1to3[new_aa]
            except KeyError:
                pass  # Keep original name if new AA is unknown
        new_chain.add(new_res)

    def _flush_insertion(self, block, anchor_prev_idx, anchor_next_idx):
        """Handle a completed insertion block."""
        if not block:
            return []
        n_ins = len(block)
        prev_ca = (
            self.ca_positions[anchor_prev_idx] if anchor_prev_idx is not None else None
        )
        next_ca = (
            self.ca_positions[anchor_next_idx] if anchor_next_idx is not None else None
        )

        start_ca = (
            prev_ca
            if prev_ca is not None
            else (
                next_ca - np.array([n_ins * 3.8, 0, 0])
                if next_ca is not None
                else np.array([0.0, 0.0, 0.0])
            )
        )
        frag = build_peptide_fragment(n_ins, start_ca=start_ca)
        anchor_next_coords = (
            {"N": None, "CA": next_ca, "C": None}
            if next_ca is not None
            else {"CA": None}
        )

        if prev_ca is not None and next_ca is not None:
            ccd_close_loop(frag, anchor_next_coords)

        bio_res = residues_to_biopy(frag, start_resseq=1, resname="GLY")
        for i, aa in enumerate(block):
            try:
                three = protein_letters_1to3[aa]
            except KeyError:
                three = "GLY"
            bio_res[i].resname = three
        return bio_res

    def build_modified_structure(self, aligned_ref, aligned_new):  # pylint: disable=too-many-locals
        """Build a new structure based on the alignment of the reference and new sequences."""
        model = self.structure[0]
        if self.chain_id not in model:
            raise KeyError(f"Chain {self.chain_id} not found")
        orig_chain = model[self.chain_id]
        _, ref_residues = get_chain_sequence(orig_chain)

        new_model = Model.Model(0)
        new_chain = Chain.Chain(self.chain_id)

        ref_index = 0
        new_res_id = 1
        self.ca_positions = [
            np.array(r["CA"].get_vector().get_array(), dtype=float)
            if "CA" in r
            else None
            for r in ref_residues
        ]
        insertion_block = []

        for a_ref, a_new in zip(aligned_ref, aligned_new):
            if a_ref != "-" and a_new != "-":
                if insertion_block:
                    new_res_list = self._flush_insertion(
                        insertion_block,
                        ref_index - 1 if ref_index > 0 else None,
                        ref_index if ref_index < len(ref_residues) else None,
                    )
                    for r in new_res_list:
                        r.id = (" ", new_res_id, " ")
                        new_chain.add(r)
                        new_res_id += 1
                    insertion_block = []

                orig = ref_residues[ref_index]
                self._clone_and_mutate_residue(new_chain, orig, new_res_id, a_new)
                ref_index += 1
                new_res_id += 1
            elif a_ref != "-" and a_new == "-":
                ref_index += 1  # Deletion
            elif a_ref == "-" and a_new != "-":
                insertion_block.append(a_new)  # Insertion

        if insertion_block:
            new_res_list = self._flush_insertion(
                insertion_block, ref_index - 1 if ref_index > 0 else None, None
            )
            for r in new_res_list:
                r.id = (" ", new_res_id, " ")
                new_chain.add(r)
                new_res_id += 1

        new_model.add(new_chain)
        new_structure = Structure.Structure("modified")
        new_structure.add(new_model)
        return new_structure

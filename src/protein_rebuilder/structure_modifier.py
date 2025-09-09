"""
Modify a protein structure based on a new sequence.
"""
from Bio.Data.PDBData import protein_letters_3to1, protein_letters_1to3
from Bio.PDB import Model, Chain, Residue, Atom
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
            aa = 'X'
        seq.append(aa)
        res_list.append(res)
    return ''.join(seq), res_list

class StructureModifier:
    """
    Modifies a protein structure based on a new sequence.
    """
    def __init__(self, structure, chain_id='A'):
        self.structure = structure
        self.chain_id = chain_id

    def build_modified_structure(self, aligned_ref, aligned_new):
        """
        Build a new structure based on the alignment of the reference and new sequences.
        """
        model = self.structure[0]
        if self.chain_id not in model:
            raise KeyError(f"Chain {self.chain_id} not found")
        orig_chain = model[self.chain_id]
        ref_seq, ref_residues = get_chain_sequence(orig_chain)

        new_model = Model.Model(0)
        new_chain = Chain.Chain(self.chain_id)

        ref_index = 0
        new_res_id = 1
        # CA positions cache
        ca_positions = []
        for r in ref_residues:
            ca = None
            if 'CA' in r:
                ca = np.array(r['CA'].get_vector().get_array(), dtype=float)
            ca_positions.append(ca)

        # We'll collect inserted residue blocks to perform CCD closure when we reach next anchor
        insertion_block = []  # list of single-letter aa for current insertion

        def flush_insertion(block, anchor_prev_idx, anchor_next_idx):
            # called when we have a completed insertion block (one or more residues)
            if not block:
                return []
            n_ins = len(block)
            # anchors: previous CA if available, next CA if available
            prev_ca = ca_positions[anchor_prev_idx] if anchor_prev_idx is not None else None
            next_ca = ca_positions[anchor_next_idx] if anchor_next_idx is not None else None

            # Build initial fragment
            start_ca = prev_ca if prev_ca is not None else (next_ca - np.array([n_ins*3.8,0,0]) if next_ca is not None else np.array([0.,0.,0.]))
            frag = build_peptide_fragment(n_ins, start_CA=start_ca)
            # Prepare anchor coords as dictionaries
            anchor_prev_coords = {'CA': prev_ca} if prev_ca is not None else None
            anchor_next_coords = {'N': None, 'CA': next_ca, 'C': None} if next_ca is not None else {'CA': None}

            # If both anchors are present, run CCD to try to close
            if prev_ca is not None and next_ca is not None:
                ccd_close_loop(frag, anchor_prev_coords, anchor_next_coords)
                # we ignore ok/fail for now; proceed with built fragment
            # Convert to biopy residues
            bio_res = residues_to_biopy(frag, start_resseq=1, resname='GLY')
            # rename residues according to block amino acids (change resname)
            for i, aa in enumerate(block):
                try:
                    three = protein_letters_1to3[aa]
                except KeyError:
                    three = 'GLY'
                bio_res[i].resname = three
            return bio_res

        # iterate alignment
        for a_ref, a_new in zip(aligned_ref, aligned_new):
            if a_ref != '-' and a_new != '-':
                # substitution or identical
                # flush any pending insertion
                if insertion_block:
                    # previous anchor is ref_index-1, next anchor is ref_index
                    new_res_list = flush_insertion(insertion_block, ref_index-1 if ref_index>0 else None, ref_index if ref_index < len(ref_residues) else None)
                    for r in new_res_list:
                        # set proper resid numbering
                        rid = (' ', new_res_id, ' ')
                        r.id = rid
                        new_chain.add(r)
                        new_res_id += 1
                    insertion_block = []

                orig = ref_residues[ref_index]
                # clone atoms
                rid = (' ', new_res_id, ' ')
                new_res = Residue.Residue(rid, orig.get_resname(), orig.get_segid())
                for atom in orig:
                    at = Atom.Atom(atom.get_name(), atom.get_vector().get_array(), atom.get_bfactor(),
                                   atom.get_occupancy(), atom.get_altloc(), atom.get_fullname(),
                                   atom.get_serial_number(), element=atom.element)
                    new_res.add(at)
                # rename if needed
                if a_ref != a_new:
                    try:
                        new_res.resname = protein_letters_1to3[a_new]
                    except KeyError:
                        pass
                new_chain.add(new_res)
                ref_index += 1
                new_res_id += 1

            elif a_ref != '-' and a_new == '-':
                # deletion: skip ref residue
                # flush insertion if any (nothing to flush in deletion)
                ref_index += 1
            elif a_ref == '-' and a_new != '-':
                # insertion relative to reference: collect AA and continue
                insertion_block.append(a_new)
            else:
                # gap-gap (shouldn't happen)
                continue

        # flush tail insertion if any (no next anchor)
        if insertion_block:
            new_res_list = flush_insertion(insertion_block, ref_index-1 if ref_index>0 else None, None)
            for r in new_res_list:
                rid = (' ', new_res_id, ' ')
                r.id = rid
                new_chain.add(r)
                new_res_id += 1
            insertion_block = []

        new_model.add(new_chain)
        from Bio.PDB import Structure
        new_structure = Structure.Structure('modified')
        new_structure.add(new_model)
        return new_structure
# protein_rebuilder/structure_modifier.py
from Bio.PDB import Model, Chain, Residue, Atom
from Bio.PDB.Polypeptide import three_to_one, one_to_three
import numpy as np

# Helper to get sequence and mapping from structure
def get_chain_sequence(chain):
    seq = []
    res_list = []
    for res in chain:
        hetfield, resseq, icode = res.get_id()
        # skip HETATM that are not residues with standard aa names
        resname = res.get_resname()
        try:
            aa = three_to_one(resname)
        except Exception:
            aa = 'X'
        seq.append(aa)
        res_list.append(res)
    return ''.join(seq), res_list

class StructureModifier:
    """
    Map aligned sequences onto the reference chain, perform:
    - substitutions: rename residue
    - deletions: remove residue
    - insertions: create placeholder residue with backbone atoms interpolated
    """

    def __init__(self, structure, chain_id='A'):
        """
        structure: Bio.PDB structure
        chain_id: which chain to operate on
        """
        self.structure = structure
        self.chain_id = chain_id

    def build_modified_structure(self, aligned_ref, aligned_new):
        """
        aligned_ref: aligned reference sequence (with '-')
        aligned_new: aligned new sequence (with '-')
        Returns a new Bio.PDB structure (shallow copy) with residues modified.
        """
        # Find original chain
        model = self.structure[0]
        if self.chain_id not in model:
            raise KeyError(f"Chain {self.chain_id} not found in model")
        orig_chain = model[self.chain_id]

        ref_seq, ref_residues = get_chain_sequence(orig_chain)
        # We'll create a new model/chain
        new_model = Model.Model(0)
        new_chain = Chain.Chain(self.chain_id)

        # iterate through aligned sequences and map ref residues accordingly
        ref_index = 0  # index into ref_residues
        new_res_id = 1

        # Last kept coordinates for interpolation for insertions
        last_ca_pos = None
        next_ca_pos_for_interp = None

        # Precompute CA positions list for interpolation if needed
        ca_positions = [None] * len(ref_residues)
        for i, r in enumerate(ref_residues):
            if 'CA' in r:
                ca_positions[i] = np.array(r['CA'].get_vector(), dtype=float)

        # Helper to create a new placeholder residue with N, CA, C, O or only CA (preferable)
        def make_placeholder_residue(resname_three, pos_ca, resseq):
            # create a new residue object
            res_id = (' ', resseq, ' ')
            new_res = Residue.Residue(res_id, resname_three, '')
            # we will add a CA atom at least
            ca = Atom.Atom('CA', pos_ca, 0.0, 1.0, ' ', 'CA', resseq, element='C')
            new_res.add(ca)
            # optionally one could add N, C, O positions by offsetting along a small vector.
            return new_res

        # iterate positions of aligned sequences
        for a_ref_char, a_new_char in zip(aligned_ref, aligned_new):
            if a_ref_char != '-' and a_new_char != '-':
                # substitution or identical -> take ref residue coords, rename if needed
                orig_res = ref_residues[ref_index]
                # clone residue: create a new Residue and copy atoms
                res_id = (' ', new_res_id, ' ')
                new_res = Residue.Residue(res_id, orig_res.get_resname(), orig_res.get_segid())
                # copy atoms (positions preserved)
                for atom in orig_res:
                    # copy atom
                    at = Atom.Atom(atom.get_name(), atom.get_vector().get_array(), atom.get_bfactor(),
                                   atom.get_occupancy(), atom.get_altloc(), atom.get_fullname(),
                                   atom.get_serial_number(), element=atom.element)
                    new_res.add(at)
                # rename residue if amino acid identity changed
                if a_ref_char != a_new_char:
                    # try to map to three letter code
                    try:
                        new_three = one_to_three(a_new_char)
                        new_res.resname = new_three
                    except Exception:
                        # keep original name if unknown
                        pass
                new_chain.add(new_res)
                # move forward
                last_ca_pos = ca_positions[ref_index]
                ref_index += 1
                new_res_id += 1

            elif a_ref_char != '-' and a_new_char == '-':
                # deletion in new sequence -> skip this residue (do not add)
                ref_index += 1
                # update last_ca_pos
                last_ca_pos = ca_positions[ref_index-1] if (ref_index-1) < len(ca_positions) else last_ca_pos

            elif a_ref_char == '-' and a_new_char != '-':
                # insertion in new sequence relative to reference -> create a placeholder residue
                # Interpolate CA between last_ca_pos and next CA (if available)
                # find next non-gap CA
                next_index = ref_index
                next_ca = None
                while next_index < len(ca_positions) and ca_positions[next_index] is None:
                    next_index += 1
                if next_index < len(ca_positions):
                    next_ca = ca_positions[next_index]
                if last_ca_pos is None and next_ca is None:
                    # fallback: put at origin
                    pos_ca = np.array([0.0, 0.0, 0.0])
                elif last_ca_pos is None:
                    pos_ca = next_ca - np.array([1.5, 0.0, 0.0])
                elif next_ca is None:
                    pos_ca = last_ca_pos + np.array([1.5, 0.0, 0.0])
                else:
                    # place between last_ca_pos and next_ca proportionally
                    pos_ca = (last_ca_pos + next_ca) / 2.0

                # create residue with appropriate three-letter code
                try:
                    three = one_to_three(a_new_char)
                except Exception:
                    three = 'UNK'
                placeholder = make_placeholder_residue(three, pos_ca, new_res_id)
                new_chain.add(placeholder)
                last_ca_pos = pos_ca
                new_res_id += 1
                # ref_index remains same (since insertion relative to ref)
            else:
                # gap-gap (should not happen)
                continue

        new_model.add(new_chain)
        # build a new structure wrapper
        from Bio.PDB import Structure
        new_structure = Structure.Structure('modified')
        new_structure.add(new_model)
        return new_structure

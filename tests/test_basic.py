import os
from protein_rebuilder.io_helpers import read_pdb_file
from protein_rebuilder.sequence_aligner import SequenceAligner
from protein_rebuilder.structure_modifier import StructureModifier, get_chain_sequence

def test_read_pdb():
    here = os.path.dirname(__file__)
    # use small PDB snippet shipped with tests or user should provide path via env
    # For test purposes, we check code paths only when a PDB is present
    test_pdb = os.environ.get("PR_TEST_PDB")
    if not test_pdb:
        return
    struct = read_pdb_file(test_pdb)
    assert struct is not None

def test_alignment_simple():
    ref = "MKT"
    new = "MKT"
    a_ref,a_new = SequenceAligner(ref,new).align()
    assert a_ref.replace('-','') == ref

def test_structure_modifier_smoke():
    test_pdb = os.environ.get("PR_TEST_PDB")
    if not test_pdb:
        return
    struct = read_pdb_file(test_pdb)
    chain = struct[0]['A']
    ref_seq, _ = get_chain_sequence(chain)
    # create new seq with a small insertion
    new_seq = ref_seq[:10] + "GG" + ref_seq[10:]
    a_ref,a_new = SequenceAligner(ref_seq, new_seq).align()
    mod = StructureModifier(struct, chain_id='A')
    new_struct = mod.build_modified_structure(a_ref,a_new)
    assert new_struct is not None
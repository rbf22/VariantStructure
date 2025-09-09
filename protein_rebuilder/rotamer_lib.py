"""
Simple backbone-independent rotamer library.

For each residue type we list a few common chi-angle tuples (in degrees).
This is small: rotamers are multiples of 60 degrees (a coarse Dunbrack-like set).
You can replace or expand this by loading a full Dunbrack library.

Format:
ROTAMERS['ARG'] = [ (chi1, chi2, chi3, chi4), ... ]
"""

ROTAMERS = {
    'ALA': [ () ],
    'VAL': [ (60,), (180,), (-60,) ],
    'LEU': [ (60,60), (180,60), (-60,180) ],
    'ILE': [ (60,60), (180,60) ],
    'SER': [ (60,), (180,), (-60,) ],
    'THR': [ (60,), (180,), (-60,) ],
    'PHE': [ (60, 180), (180, -60), (-60, 60) ],
    'TYR': [ (60, 180), (180, -60), (-60, 60) ],
    'TRP': [ (60, 180), (180, -60) ],
    'ASP': [ (60,180), (180,-60) ],
    'GLU': [ (60,180), (180,-60), (-60,60) ],
    'ASN': [ (60, 180), (180, -60) ],
    'GLN': [ (60, 180), (180, -60), (-60, 60) ],
    'HIS': [ (60, 180), (180, -60) ],
    'LYS': [ (60, 60, 60), (180, 60, -60), (-60, -60, 60) ],
    'ARG': [ (60, 60, 60, 60), (180, 60, 60, -60) ],
    'PRO': [ ( -60,), ( -30,) ],
    'CYS': [ (60,), (180,), (-60,) ],
    'GLY': [ () ]
}

def get_rotamers_for_resname(resname_three):
    import sys
    from Bio.PDB.Polypeptide import three_to_one
    try:
        aa1 = three_to_one(resname_three)
    except Exception:
        aa1 = 'X'
    # map one-letter to three-letter keys used above
    # convert one-letter to three-letter standard
    from Bio.PDB.Polypeptide import one_to_three
    try:
        three = one_to_three(aa1)
    except Exception:
        three = resname_three
    return ROTAMERS.get(three, [()])
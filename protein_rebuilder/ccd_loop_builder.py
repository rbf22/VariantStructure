"""
CCD-based loop builder:

Given two anchor residues (prev_res, next_res) in 3D and a number of
inserted residues n_ins, create an initial backbone for the inserted
residues using ideal peptide geometry and then apply CCD to adjust
phi angles of the inserted residues so that the C-terminal anchor
(next_res N atom) is reached (loop closed).

The approach:
- Build an initial chain of N/CA/C/O atoms for each inserted residue using ideal bond lengths/angles.
- Use CCD altering phi dihedrals of inserted residues to move the chain tip to target.
- Return list of Residue objects (Bio.PDB Residue) with N, CA, C, O atoms set.

This implementation is intentionally compact and focuses on correctness for small loops.
"""

from Bio.PDB import Residue, Atom
import numpy as np
from math import cos, sin, radians, acos
from Bio.PDB.Polypeptide import one_to_three

# Standard backbone geometry (angstroms/degrees)
BOND_N_CA = 1.458  # N-CA
BOND_CA_C = 1.525  # CA-C
BOND_C_N = 1.329   # C-N (peptide)
BOND_C_O = 1.229   # C-O
ANGLE_C_N_CA = 121.7  # degrees (C-N-CA)
ANGLE_N_CA_C = 110.4  # degrees
ANGLE_CA_C_N = 116.6  # degrees
OMEGA = 180.0  # trans peptide

# Utility vector math
def norm(v): return v / np.linalg.norm(v)
def length(v): return np.linalg.norm(v)
def dihedral(p0, p1, p2, p3):
    # returns dihedral angle in degrees between 4 points
    b0 = -1.0 * (p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    # normalize b1
    b1 /= np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))

def set_dihedral(p0, p1, p2, p3, new_angle_deg):
    # rotate p3 around axis p1-p2 so that dihedral(p0,p1,p2,p3) == new_angle_deg
    current = dihedral(p0, p1, p2, p3)
    delta = np.radians(new_angle_deg - current)
    axis = p2 - p1
    axis = axis / np.linalg.norm(axis)
    # Rodrigues rotation of p3 around axis passing through p2
    v = p3 - p2
    ca = np.cos(delta)
    sa = np.sin(delta)
    v_rot = v * ca + np.cross(axis, v) * sa + axis * np.dot(axis, v) * (1 - ca)
    return p2 + v_rot

def build_peptide_fragment(n_residues, start_CA=np.array([0.0,0.0,0.0]), direction=np.array([1.5,0.0,0.0])):
    """
    Create an initial straight fragment of n_residues backbone coordinates (N,CA,C,O) using ideal geometry.
    Returns list of residues where each residue is dict of atom name -> 3-vector.
    """
    residues = []
    # We'll place CA atoms along direction spaced by ~3.8 A
    ca_sep = 3.8
    for i in range(n_residues):
        ca_pos = start_CA + direction * (i * ca_sep)
        # approximate N and C positions relative to CA along the chain direction
        n_pos = ca_pos - np.array([1.45, 0.2, 0.0])
        c_pos = ca_pos + np.array([1.52, -0.2, 0.0])
        o_pos = c_pos + np.array([0.5, -0.4, 0.0])
        residues.append({'N': n_pos, 'CA': ca_pos, 'C': c_pos, 'O': o_pos})
    return residues

def residues_to_biopy(res_dicts, start_resseq=1, resname='GLY'):
    """
    Convert residue dictionaries to Bio.PDB Residue objects (with Atom objects).
    Only backbone atoms are added; element names are inferred.
    """
    from Bio.PDB import Residue, Atom
    bio_residues = []
    for i, rd in enumerate(res_dicts):
        rid = (' ', start_resseq + i, ' ')
        res = Residue.Residue(rid, resname, '')
        for name in ['N','CA','C','O']:
            pos = rd[name]
            atom = Atom.Atom(name, pos, 0.0, 1.0, ' ', name, i+1, element=name[0])
            res.add(atom)
        bio_residues.append(res)
    return bio_residues

def ccd_close_loop(insert_coords, anchor_prev_coords, anchor_next_coords, max_iter=200, tol=1e-3):
    """
    insert_coords: list of residue dicts for inserted residues (with 'N','CA','C','O' arrays)
    anchor_prev_coords: dictionary of anchor residue atoms (the residue before insertion)
    anchor_next_coords: dictionary for the residue after insertion

    We attempt to adjust phi dihedrals (rotation around C(i-1)-N(i)-CA(i)-C(i)) of the inserted residues
    via CCD to bring the last built N atom close to anchor_next_coords['N'] (or next anchor N).
    This is a simplified CCD implementation operating on points in space.
    """
    # Build points list: we will treat residues as chain of CA positions for CCD control
    # For better control, we rotate about virtual bonds affecting positions after that bond.
    # We'll operate by adjusting phi (rotation of N-CA-C atoms) — simplified:
    target = anchor_next_coords['N']
    for it in range(max_iter):
        end_CA = insert_coords[-1]['CA']
        dist = length(end_CA - target)
        if dist < tol:
            return True
        # CCD: iterate from last movable joint backwards
        for j in range(len(insert_coords)-1, -1, -1):
            # joint position: CA_j (we choose CA as pivot)
            pivot = insert_coords[j]['CA']
            end = insert_coords[-1]['CA']
            # vector from pivot to end and to target
            v_end = end - pivot
            v_tgt = target - pivot
            if length(v_end) < 1e-6 or length(v_tgt) < 1e-6:
                continue
            v_end_u = v_end / length(v_end)
            v_tgt_u = v_tgt / length(v_tgt)
            # compute rotation axis and angle
            axis = np.cross(v_end_u, v_tgt_u)
            axis_norm = length(axis)
            if axis_norm < 1e-8:
                continue
            axis = axis / axis_norm
            cosang = np.dot(v_end_u, v_tgt_u)
            cosang = max(-1.0, min(1.0, cosang))
            angle = np.arccos(cosang)
            # rotate all points after joint around axis by angle (Rodrigues)
            for k in range(j+1, len(insert_coords)):
                for a in ['N','CA','C','O']:
                    p = insert_coords[k][a] - pivot
                    ca = np.cos(angle)
                    sa = np.sin(angle)
                    p_rot = p * ca + np.cross(axis, p) * sa + axis * np.dot(axis, p) * (1 - ca)
                    insert_coords[k][a] = pivot + p_rot
        # next iteration
    return False
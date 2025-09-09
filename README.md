# VariantStructure
A small, self-contained Python package that uses PDBFixer to rebuild missing atoms/sidechains and OpenMM to minimize the structure for a new sequence and a reference structure.

 - Aligns a new sequence to the sequence present in a reference PDB (via Biopython pairwise2),
 - Maps mutations/insertions/deletions onto the reference chain,
 - Creates placeholder backbone atoms for inserted residues (simple linear interpolation of backbone positions),
 - Uses PDBFixer to build missing atoms/sidechains and add hydrogens,
 - Uses OpenMM to perform restrained minimization to relieve clashes.

We provide an H-Ras example PDB for the demo (PDB ID 4EFL, a GppNHp-bound H-Ras structure). (Reference: RCSB PDB entry for H-Ras; example PDB: 4EFL). 

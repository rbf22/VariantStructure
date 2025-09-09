# Protein Rebuilder

A small, self-contained Python package that uses PDBFixer to rebuild missing atoms/sidechains and OpenMM to minimize the structure for a new sequence and a reference structure.

- Aligns a new sequence to the sequence present in a reference PDB (via Biopython pairwise2),
- Maps mutations/insertions/deletions onto the reference chain,
- Creates placeholder backbone atoms for inserted residues (simple linear interpolation of backbone positions),
- Uses PDBFixer to build missing atoms/sidechains and add hydrogens,
- Uses OpenMM to perform restrained minimization to relieve clashes.

## Installation

To install the package, clone the repository and install it using pip:

```bash
git clone https://github.com/example/protein-rebuilder.git
cd protein-rebuilder
pip install .
```

## Usage

The package provides a command-line tool `protein-rebuilder` to rebuild a protein structure.

```bash
protein-rebuilder --pdb <reference_pdb_file> --new-seq <new_sequence_file> --out <output_pdb_file>
```

### Arguments

- `--pdb`: Path to the reference PDB file.
- `--new-seq`: Path to the new sequence file (in FASTA or plain text format).
- `--out`: Path to the output PDB file.
- `--chain`: The chain to operate on (default: A).
- `--min-iterations`: Number of minimization iterations (default: 500).
- `--mc-iterations`: Number of rotamer Monte Carlo iterations (default: 200).

### Example

An example of how to use the tool is provided in `example/example-hras.py`. To run the example, you first need to install the development dependencies:

```bash
pip install -e .[dev]
```

Then you can run the example:

```bash
python example/example-hras.py
```

## Running tests

To run the tests, first install the development dependencies:

```bash
pip install -e .[dev]
```

Then run pytest:

```bash
pytest
```

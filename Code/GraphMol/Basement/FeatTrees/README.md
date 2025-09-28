# Feature Tree Utilities

This directory collects experimental tooling for constructing and comparing
feature trees derived from RDKit molecules.  The utilities operate on the
Boost-based feature-tree graph produced by `FeatTrees::molToBaseTree()` and
provide both command-line and Python access to fragment-level descriptors and
similarity scoring.

## Contents

- `FeatTree.h` / `FeatTree.cpp`: definitions of the base feature-tree graph and
  helpers that cluster fused rings and connector atoms.
- `FeatTreeUtils.h` / `FeatTreeUtils.cpp`: helper routines used while building
  the base tree.
- `FeatTreeSimilarity.h` / `FeatTreeSimilarity.cpp`: fragment profiling and
  similarity calculation utilities shared by the command-line program and the
  Python bindings.
- `ftrees_full.cpp`: a standalone command-line tool that evaluates pairwise
  feature-tree similarities for SMILES inputs.
- `testFeatTrees.cpp`: legacy exploratory tests for the feature-tree
  construction logic.

## Command-line usage

The `ftrees_full` executable expects two or more SMILES strings.  Each SMILES
is parsed, sanitized, converted into a feature tree, and compared to every other
input molecule using the fragment-level similarity heuristic.

```bash
ftrees_full "c1ccccc1O" "CCN(CC)CC"
```

Each pair of molecules is reported with a similarity score in the range `[0, 1]`.

## Python bindings

The `rdftrees` Python extension exposes the same profiling and scoring logic.
It can be imported from `rdkit.Chem` once RDKit is built with Python bindings:

```python
from rdkit import Chem
from rdkit.Chem import rdftrees

mol1 = Chem.MolFromSmiles("c1ccccc1O")
mol2 = Chem.MolFromSmiles("CCN(CC)CC")

profiles = rdftrees.build_fragment_profiles(mol1)
score = rdftrees.compare_molecules(mol1, mol2)
```

- `build_fragment_profiles(mol, sanitize=True)` returns a list of dictionaries,
  one per fragment, summarizing the coarse physicochemical properties of the
  fragment.
- `compare_molecules(mol1, mol2, sanitize=True)` computes a symmetric
  similarity score by matching the best fragment pairings between both
  molecules.

All functions sanitize a working copy of the provided molecules by default.
Disable sanitization only if the molecules were already processed.

## Implementation notes

The fragment descriptors capture the number of heavy atoms (volume), ring and
aromatic membership flags, counts of hydrogen-bond donors and acceptors, a
simple amide-like pattern, and the fraction of carbon atoms as a proxy for
hydrophobicity.  Similarities are computed by comparing these descriptors with a
set of heuristic component scores and averaging the best-matching fragment
pairs between the molecules.


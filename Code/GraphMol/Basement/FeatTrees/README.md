# Feature Trees

Modern feature-tree construction and similarity evaluation tools for RDKit. The
implementation builds deterministic graph abstractions that capture ring
systems, connector motifs, terminal branches, and optional pharmacophoric
feature groups.  Each node stores sorted atom indices alongside derived summary
statistics which can be fed into downstream scoring algorithms.

## Contents

- `FeatTree.h` / `FeatTree.cpp`: public API and implementation of the base-tree
  builder, the feature-tree transformer, canonicalisation and JSON utilities,
  together with the weighted Jaccard similarity calculation.
- `FeatTreeUtils.h` / `FeatTreeUtils.cpp`: factored helper routines used during
  the base-tree construction pass.  These remain available for unit testing even
  though `molToBaseTree()` now performs most steps internally.
- `FeatTreeSimilarity.h` / `FeatTreeSimilarity.cpp`: tree edit-distance
  approximation helpers that piggy-back on the similarity scores.
- `FeatTrees/Wrap/rdFeatTrees.cpp`: Python bindings exposing the new
  functionality as `rdkit.Chem.rdFeatTrees`.
- `testFeatTreeConstruction.cpp`, `testFeatTreeSimilarity.cpp`: unit tests that
  validate construction invariants, determinism, feature annotations, and the
  similarity interface.
- `ftrees_full.cpp`: legacy command-line driver retained for manual regression
  checks.

## Quick start

```cpp
#include <GraphMol/Basement/FeatTrees/FeatTree.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

using namespace RDKit;
using namespace RDKit::FeatTrees;

ROMol mol(*SmilesToMol("c1ccc(CO)cc1"));
FeatTreeParams params;
params.annotateFeatures = true;
auto tree = molToFeatTree(mol, params);
std::cout << featTreeToJSON(*tree) << std::endl;
```

Python bindings expose a similar workflow:

```python
from rdkit import Chem
from rdkit.Chem import rdFeatTrees

mol = Chem.MolFromSmiles("c1ccc(O)cc1")
params = rdFeatTrees.FeatTreeParams(annotateFeatures=True)
tree = rdFeatTrees.MolToFeatTree(mol, params)
print(tree.ToJSON())
```

## Determinism and testing

Both the C++ and Python builders canonicalise node ordering before returning a
feature tree.  The regression tests in this directory construct multiple
reference trees, exercise parameter toggles (path compression, branch grouping,
feature annotation), and ensure that serialised JSON descriptions remain stable
across repeated runs.  Similarity tests guard the weighted Jaccard scores and
edit-distance approximation helpers.

## Similarity scoring

`calcFeatTreeSimilarity()` implements a weighted Jaccard index over node
signatures.  Each signature summarises the node kind, atom counts, aromatic and
hetero-atom statistics, and selected flags (donor/acceptor).  The same scoring
machinery feeds the edit-distance approximation in `FeatTreeSimilarity.cpp`.

## Feature annotations

When `FeatTreeParams::annotateFeatures` is enabled the transformer attempts to
load the default RDKit feature definition file (`BaseFeatures.fdef`).  If the
file is unavailable the code falls back to a small heuristic detector that tags
hetero atoms.  Feature-group nodes coexist with structural nodes and are marked
with the `FeatTreeNodeKind::FeatureGroup` kind for easy filtering.

## Python module

The `rdkit.Chem.rdFeatTrees` module exposes:

- `FeatTreeParams`: keyword-configurable parameter object with a helpful `repr`.
- `MolToFeatTree(mol, params=None, as_base_tree=False)`: builds either the base
  tree or the fully transformed tree.
- `BaseTreeToFeatTree(base_tree, params=None, mol)`: applies the transformation
  stages to an existing base tree instance.
- `CalcFeatTreeSimilarity(...)`: overloaded for molecule or `FeatTree` inputs.
- `FeatTree` class: `.GetNodes()`, `.GetEdges()`, `.ToJSON()`, and a concise
  `repr` string.

All Python node dictionaries mirror the C++ `FeatTreeNodeData` structure so that
flags and summary statistics stay accessible from scripts.

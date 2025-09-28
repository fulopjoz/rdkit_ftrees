Feature Trees
=============

The feature-tree machinery builds concise graph abstractions over molecules so
that ring systems, connector motifs, and optional pharmacophoric features can be
reasoned about using graph algorithms.  The API is implemented in
``GraphMol/Basement/FeatTrees`` and is available in both C++ and Python.

Overview
--------

A feature tree is an undirected graph with typed nodes:

* ``RingSystem`` and ``FusedRingSystem`` nodes collect fused cycles.
* ``Connector`` nodes represent non-ring atoms that bridge structural groups.
* ``BranchGroup`` nodes aggregate terminal branches when requested.
* ``FeatureGroup`` nodes (optional) capture pharmacophoric features.
* ``ZeroNode`` nodes (optional) break cycles introduced by non-ring atoms.

Each node stores the sorted list of atom indices that it represents alongside
summary statistics (aromatic count, hetero count, donor/acceptor flags, and ring
size ranges).  Edges record whether either endpoint belongs to a ring-based
node.

The resulting structure is deterministic.  Canonicalisation sorts nodes using a
priority key so that serialisation and comparisons remain reproducible.

Simple Example
--------------

::

        Ring(Fused)      Connector      BranchGroup
             o---------------o---------------o
             |                               |
             o-------------------------------o

C++ usage::

    using namespace RDKit;
    using namespace RDKit::FeatTrees;

    ROMol mol(*SmilesToMol("c1ccc(OCCN)cc1"));
    FeatTreeParams params;
    params.annotateFeatures = true;
    params.maxBranchGroupSize = 2;
    auto tree = molToFeatTree(mol, params);
    std::cout << featTreeToJSON(*tree) << std::endl;

Python usage::

    from rdkit import Chem
    from rdkit.Chem import rdFeatTrees

    mol = Chem.MolFromSmiles("c1ccc(OCCN)cc1")
    params = rdFeatTrees.FeatTreeParams(annotateFeatures=True,
                                        maxBranchGroupSize=2)
    tree = rdFeatTrees.MolToFeatTree(mol, params)
    print(tree.ToJSON())

Parameters
----------

``FeatTreeParams`` controls the construction and post-processing stages:

=====================  ================================================
Parameter              Meaning
=====================  ================================================
``compressPaths``      Collapse chains of degree-2 connectors into a single
                       node (default ``True``).
``mergeBranchGroups``  Convert terminal connectors with small atom counts into
                       ``BranchGroup`` nodes (default ``True``).
``maxBranchGroupSize`` Maximum number of atoms to aggregate in a branch group
                       (default ``3``).
``annotateFeatures``   Add pharmacophoric ``FeatureGroup`` nodes.  The code
                       attempts to load ``BaseFeatures.fdef`` and falls back to
                       a simple hetero-atom heuristic if unavailable.
``includeZeroNodes``   Insert ``ZeroNode`` breakpoints along non-ring cycles
                       (default ``False`` because the new edge flags normally
                       suffice).
``canonicalize``       Sort nodes deterministically before returning the graph.
``ringWeight``         Weight multiplier used by the weighted Jaccard similarity
                       for ring-based nodes.
``connectorWeight``    Similarity weight for connector and branch nodes.
``featureGroupWeight`` Similarity weight for feature nodes.
=====================  ================================================

Similarity
----------

``calcFeatTreeSimilarity`` computes a weighted Jaccard index over node
signatures.  Each signature contains the node kind, atom count, aromatic and
hetero counts, ring size range, and summarised flags.  The same values feed the
``calcFeatTreeEditDistanceApprox`` helper, which returns ``(1 - similarity)
* average_node_count`` as a fast tree-edit-distance surrogate.

Performance Notes
-----------------

* Ring information is computed once per molecule using RDKit's ``RingInfo``
  caches.
* Atom vectors are stored as sorted ``std::vector<unsigned int>`` instances.
  ``finalizeAtomVector`` enforces determinism.
* Path compression and canonicalisation operate in ``O(V log V + E)`` time.
* JSON serialisation is intended for diagnostics and regression tests, not for
  high-volume interchange.

The unit tests ``graphmolFeatTreeConstruction`` and
``graphmolFeatTreeSimilarity`` cover the invariants, determinism checks, feature
annotation toggles, and similarity bounds.

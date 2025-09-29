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
    validateFeatTree(*tree, params, /*allowZeroNodes=*/true);
    std::cout << featTreeToJSON(*tree, params) << std::endl;

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
``similarityAutoThreshold`` Node-count threshold used when
                       ``FeatTreeSimilarityMethod::Auto`` selects between
                       weighted Jaccard and the edit-distance approximation.
=====================  ================================================

Similarity
----------

``calcFeatTreeSimilarity`` computes a weighted Jaccard index over node
signatures.  Each signature contains the node kind, atom count, aromatic and
hetero counts, ring size range, and summarised flags.  The similarity score is
``sum(min(w_i^A, w_i^B)) / sum(max(w_i^A, w_i^B))`` over the matching signature
weights ``w_i``.  The same values feed the ``calcFeatTreeEditDistanceApprox``
helper, which returns ``(1 - similarity) * average_node_count`` as a fast
tree-edit-distance surrogate.  The ``FeatTreeSimilarityMethod`` enum controls
which algorithm is applied; the ``Auto`` mode uses the edit-distance surrogate
for small graphs and falls back to weighted Jaccard otherwise.  Both algorithms
gracefully handle empty graphs, returning a perfect score when both inputs are
empty.

Validation and hashing
----------------------

``validateFeatTree`` enforces internal invariants: node atom indices remain
sorted, non-zero nodes are non-empty, structural nodes form a disjoint partition
of atom indices, edges avoid duplicates, ring-end counts are within ``[0, 2]``
and canonical ordering is honoured.  ``validateParams`` guards the construction
parameters.  ``hashFeatTree`` produces a 64-bit canonical fingerprint that can
be used for caching or as a determinism check.

Performance Notes
-----------------

* Ring information is computed once per molecule using RDKit's ``RingInfo``
  caches.
* Atom vectors are stored as sorted ``std::vector<unsigned int>`` instances.
  ``finalizeAtomVector`` enforces determinism.
* Path compression and canonicalisation operate in ``O(V log V + E)`` time.
* JSON serialisation emits ``schema_version`` along with the influential
  parameter settings and is intended for diagnostics and regression tests, not
  for high-volume interchange.

The unit tests ``graphmolFeatTreeConstruction``,
``graphmolFeatTreeSimilarity``, ``graphmolFeatTreeInvariants``,
``graphmolFeatTreeHash`` and ``graphmolFeatTreeSimilarityMethods`` cover the
invariants, determinism checks, hashing, feature annotation toggles and
similarity bounds.

Python module additions
-----------------------

``rdkit.Chem.rdFeatTrees`` also exposes ``HashFeatTree`` along with the
``FeatTreeSimilarityMethod`` enum and the ``FEATTREE_SCHEMA_VERSION`` constant.
Calling ``FeatTree.ToJSON(params)`` returns the serialised structure with the
current schema version and the influential parameter settings.

Thread safety and profiling
---------------------------

Feature trees are immutable snapshots once constructed.  As long as client code
avoids mutating the underlying boost graph instances they may be shared safely
across threads.  Defining ``FEATTREE_PROFILE`` during compilation enables light
profiling hooks (reported via ``rdDebug``) around construction, transformation,
canonicalisation and similarity evaluation.

Multi-fragment molecules
------------------------

The builders process all fragments present in the input ``ROMol``.  Workflows
that only care about the largest component should pre-select the desired
fragment before invoking the feature-tree API.

Build toggle and API stability
------------------------------

The CMake option ``RDK_BUILD_FEATTREES`` (enabled by default in this fork)
controls whether the feature-tree sources, tests and Python bindings are
compiled.  The API is currently considered *experimental* and may evolve as the
similarity methods are extended; new helpers will maintain the existing
invariants and canonicalisation behaviour.

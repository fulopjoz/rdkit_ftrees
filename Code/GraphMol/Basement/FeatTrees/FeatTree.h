//
//  Copyright (C) 2024  RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

/*! \file FeatTree.h

    \brief Public API for building, transforming and comparing feature trees.

    The feature-tree graph is a light-weight abstraction of a molecule that
    groups rings, connectors, terminal branches and optional feature groups into
    typed nodes.  The resulting representation is deterministic and can be
    serialised or compared using weighted Jaccard similarity scores.

    \since 2024.09
*/

#ifndef RD_FEATTREE_H
#define RD_FEATTREE_H

#include <RDGeneral/export.h>
#include <RDGeneral/Invariant.h>

#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/shared_ptr.hpp>

#include <cstdint>
#include <string>
#include <vector>

namespace RDKit {
class ROMol;

namespace FeatTrees {

//! \brief schema version for feature tree serialisation.
//! Increment when the JSON encoding gains backwards-incompatible changes.
constexpr uint32_t FEATTREE_SCHEMA_VERSION = 1;

//! \enum FeatTreeNodeKind
//! \brief Enumerates the supported structural node categories.
enum class FeatTreeNodeKind : unsigned char {
  RingSystem = 0,
  FusedRingSystem,
  Connector,
  BranchGroup,
  ZeroNode,
  FeatureGroup
};

//! \enum FeatTreeSimilarityMethod
//! \brief Selects the similarity scoring algorithm to apply.
enum class FeatTreeSimilarityMethod : unsigned char {
  WeightedJaccard = 0,
  ApproxEdit,
  Auto
};

//! \brief bit flags stored on nodes for derived annotations.
struct FeatTreeNodeFlags {
  enum : unsigned char {
    Aromatic = 0x01,
    HasHetero = 0x02,
    Donor = 0x04,
    Acceptor = 0x08,
    ContainsCharge = 0x10,
    Reserved1 = 0x20,
    Reserved2 = 0x40,
    Reserved3 = 0x80
  };
};

//! \brief bit flags stored on edges.
struct FeatTreeEdgeFlags {
  enum : unsigned char {
    RingBridge = 0x01,
    PathCompressed = 0x02,
    Reserved1 = 0x04,
    Reserved2 = 0x08,
    Reserved3 = 0x10,
    Reserved4 = 0x20,
    Reserved5 = 0x40,
    Reserved6 = 0x80
  };
};

//! \brief payload stored for every feature tree node.
struct FeatTreeNodeData {
  std::vector<unsigned int> atoms;
  FeatTreeNodeKind kind = FeatTreeNodeKind::Connector;
  unsigned char flags = 0u;
  unsigned char minRingSize = 0u;
  unsigned char maxRingSize = 0u;
  unsigned char aromaticAtomCount = 0u;
  unsigned char heteroAtomCount = 0u;
};

//! \brief payload stored for every feature tree edge.
struct FeatTreeEdgeData {
  unsigned char ringEndCount = 0u;
  unsigned char flags = 0u;
};

//! \brief property tag used by boost::adjacency_list for feature tree nodes.
struct FeatTreeNode_t {
  enum { num = 1027 };
  typedef boost::vertex_property_tag kind;
};

//! \brief property tag used by boost::adjacency_list for feature tree edges.
struct FeatTreeEdge_t {
  enum { num = 1028 };
  typedef boost::edge_property_tag kind;
};

//! \brief adjacency list that represents feature trees.
//!
//! Feature-tree graphs returned from the public builders are immutable; callers
//! should treat them as read-only snapshots.  All transformations return new
//! graphs to avoid aliasing surprises.
typedef boost::adjacency_list<
    boost::vecS, boost::vecS, boost::undirectedS,
    boost::property<FeatTreeNode_t, FeatTreeNodeData>,
    boost::property<FeatTreeEdge_t, FeatTreeEdgeData>> FeatTreeGraph;

//! shared pointer to a feature tree graph.
typedef boost::shared_ptr<FeatTreeGraph> FeatTreeGraphSPtr;

//! convenient aliases for property maps
typedef boost::property_map<FeatTreeGraph, FeatTreeNode_t>::type FeatTreeNodePMap;
typedef boost::property_map<FeatTreeGraph, FeatTreeEdge_t>::type FeatTreeEdgePMap;

//! \brief parameters controlling construction and post-processing.
struct RDKIT_GRAPHMOL_EXPORT FeatTreeParams {
  bool compressPaths = true;
  bool mergeBranchGroups = true;
  unsigned int maxBranchGroupSize = 3;
  bool annotateFeatures = false;
  bool includeZeroNodes = false;
  bool canonicalize = true;
  double ringWeight = 1.0;
  double connectorWeight = 1.0;
  double featureGroupWeight = 1.2;
  unsigned int similarityAutoThreshold = 32;
};

//! helper returning true if the given node is a ring-like node.
inline bool isRingNode(const FeatTreeNodeData &node) {
  return node.kind == FeatTreeNodeKind::RingSystem ||
         node.kind == FeatTreeNodeKind::FusedRingSystem;
}

//! helper returning true if the node is marked as zero node.
inline bool isZeroNode(const FeatTreeNodeData &node) {
  return node.kind == FeatTreeNodeKind::ZeroNode;
}

//! helper to set a node flag value.
inline void setNodeFlag(FeatTreeNodeData &node, unsigned char mask,
                        bool value) {
  if (value) {
    node.flags |= mask;
  } else {
    node.flags &= static_cast<unsigned char>(~mask);
  }
}

//! helper to test a flag value on a node.
inline bool getNodeFlag(const FeatTreeNodeData &node, unsigned char mask) {
  return (node.flags & mask) != 0u;
}

//! helper to set an edge flag value.
inline void setEdgeFlag(FeatTreeEdgeData &edge, unsigned char mask,
                        bool value) {
  if (value) {
    edge.flags |= mask;
  } else {
    edge.flags &= static_cast<unsigned char>(~mask);
  }
}

//! helper to test a flag value on an edge.
inline bool getEdgeFlag(const FeatTreeEdgeData &edge, unsigned char mask) {
  return (edge.flags & mask) != 0u;
}

//! \brief finalises a node atom vector ensuring determinism.
inline void finalizeAtomVector(std::vector<unsigned int> &atoms) {
  std::sort(atoms.begin(), atoms.end());
  atoms.erase(std::unique(atoms.begin(), atoms.end()), atoms.end());
}

//! \brief builds the base feature tree for the supplied molecule.
RDKIT_GRAPHMOL_EXPORT FeatTreeGraphSPtr molToBaseTree(
    const ROMol &mol, const FeatTreeParams &params = FeatTreeParams());

//! \brief transforms a base tree into an enriched feature tree.
RDKIT_GRAPHMOL_EXPORT void baseTreeToFeatTree(
    FeatTreeGraph &baseTree, const FeatTreeParams &params = FeatTreeParams(),
    const ROMol *mol = nullptr);

//! \brief convenience that builds the full feature tree in a single call.
RDKIT_GRAPHMOL_EXPORT FeatTreeGraphSPtr molToFeatTree(
    const ROMol &mol, const FeatTreeParams &params = FeatTreeParams());

//! \brief canonicalises vertex numbering to provide deterministic serialisation.
RDKIT_GRAPHMOL_EXPORT void canonicalizeFeatTree(FeatTreeGraph &graph);

//! \brief validates feature tree structural invariants.
RDKIT_GRAPHMOL_EXPORT void validateFeatTree(const FeatTreeGraph &graph,
                                            const FeatTreeParams &params,
                                            bool allowZeroNodes = false);

//! \brief validates construction parameters.
RDKIT_GRAPHMOL_EXPORT void validateParams(const FeatTreeParams &params);

//! \brief calculates a stable hash for canonical feature trees.
RDKIT_GRAPHMOL_EXPORT uint64_t hashFeatTree(const FeatTreeGraph &graph);

//! \brief serialises a feature tree to a JSON string for debugging.
RDKIT_GRAPHMOL_EXPORT std::string featTreeToJSON(
    const FeatTreeGraph &graph,
    const FeatTreeParams &params = FeatTreeParams());

//! \brief deserialises a feature tree from JSON (not yet implemented).
RDKIT_GRAPHMOL_EXPORT FeatTreeGraphSPtr featTreeFromJSON(
    const std::string &json, bool allowExperimental = false);

//! \brief Computes similarity between two molecules using feature trees.
RDKIT_GRAPHMOL_EXPORT double calcFeatTreeSimilarity(
    const ROMol &mol1, const ROMol &mol2,
    FeatTreeSimilarityMethod method = FeatTreeSimilarityMethod::WeightedJaccard,
    const FeatTreeParams &params = FeatTreeParams());

//! \brief Backwards compatible overload selecting weighted Jaccard.
inline double calcFeatTreeSimilarity(const ROMol &mol1, const ROMol &mol2,
                                     const FeatTreeParams &params) {
  return calcFeatTreeSimilarity(mol1, mol2,
                                FeatTreeSimilarityMethod::WeightedJaccard,
                                params);
}

//! \brief Computes similarity between two feature trees.
RDKIT_GRAPHMOL_EXPORT double calcFeatTreeSimilarity(
    const FeatTreeGraph &g1, const FeatTreeGraph &g2,
    FeatTreeSimilarityMethod method =
        FeatTreeSimilarityMethod::WeightedJaccard,
    const FeatTreeParams &params = FeatTreeParams());

//! \brief Backwards compatible overload selecting weighted Jaccard.
inline double calcFeatTreeSimilarity(const FeatTreeGraph &g1,
                                     const FeatTreeGraph &g2,
                                     const FeatTreeParams &params) {
  return calcFeatTreeSimilarity(g1, g2, FeatTreeSimilarityMethod::WeightedJaccard,
                                params);
}

}  // namespace FeatTrees
}  // namespace RDKit

#endif  // RD_FEATTREE_H

//
//  Copyright (C) 2024  RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "FeatTreeUtils.h"

#include <RDGeneral/BoostStartInclude.h>
#include <boost/graph/graph_traits.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/RingInfo.h>

#include <algorithm>
#include <iterator>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <vector>

namespace RDKit {
namespace FeatTrees {
namespace {

struct UnionFind {
  explicit UnionFind(std::size_t n) : parent(n), rank(n, 0) {
    std::iota(parent.begin(), parent.end(), 0);
  }
  unsigned int find(unsigned int x) {
    if (parent[x] != x) {
      parent[x] = find(parent[x]);
    }
    return parent[x];
  }
  void unite(unsigned int a, unsigned int b) {
    a = find(a);
    b = find(b);
    if (a == b) {
      return;
    }
    if (rank[a] < rank[b]) {
      parent[a] = b;
    } else if (rank[b] < rank[a]) {
      parent[b] = a;
    } else {
      parent[b] = a;
      ++rank[a];
    }
  }
  std::vector<unsigned int> parent;
  std::vector<unsigned int> rank;
};

struct AtomClassification {
  bool inRing = false;
  unsigned int ringGroup = std::numeric_limits<unsigned int>::max();
};

void accumulateNodeStats(const ROMol &mol, const std::vector<unsigned int> &atoms,
                         FeatTreeNodeData &node) {
  node.aromaticAtomCount = 0u;
  node.heteroAtomCount = 0u;
  bool donor = false;
  bool acceptor = false;
  bool charge = false;
  for (auto idx : atoms) {
    const auto at = mol.getAtomWithIdx(idx);
    if (!at) {
      continue;
    }
    if (at->getIsAromatic()) {
      ++node.aromaticAtomCount;
    }
    if (at->getAtomicNum() != 6 && at->getAtomicNum()) {
      ++node.heteroAtomCount;
    }
    donor = donor || (at->getTotalNumHs(true) > 0 && at->getAtomicNum() != 6);
    acceptor = acceptor || (at->getTotalValence() - at->getExplicitValence() > 0);
    charge = charge || (at->getFormalCharge() != 0);
  }
  setNodeFlag(node, FeatTreeNodeFlags::Aromatic, node.aromaticAtomCount > 0u);
  setNodeFlag(node, FeatTreeNodeFlags::HasHetero, node.heteroAtomCount > 0u);
  setNodeFlag(node, FeatTreeNodeFlags::Donor, donor);
  setNodeFlag(node, FeatTreeNodeFlags::Acceptor, acceptor);
  setNodeFlag(node, FeatTreeNodeFlags::ContainsCharge, charge);
}

std::vector<AtomClassification> classifyAtoms(const ROMol &mol,
                                              std::vector<std::vector<int>> &rings,
                                              std::vector<std::vector<unsigned int>> &groups,
                                              std::vector<unsigned int> &ringSizes) {
  RingInfo *ringInfo = mol.getRingInfo();
  rings = ringInfo->atomRings();
  for (auto &ring : rings) {
    std::sort(ring.begin(), ring.end());
  }
  ringSizes.resize(rings.size());
  for (unsigned int i = 0; i < rings.size(); ++i) {
    ringSizes[i] = static_cast<unsigned int>(rings[i].size());
  }
  if (rings.empty()) {
    return std::vector<AtomClassification>(mol.getNumAtoms());
  }
  UnionFind uf(rings.size());
  for (unsigned int i = 0; i < rings.size(); ++i) {
    for (unsigned int j = i + 1; j < rings.size(); ++j) {
      std::vector<int> tmp;
      std::set_intersection(rings[i].begin(), rings[i].end(), rings[j].begin(),
                            rings[j].end(), std::back_inserter(tmp));
      if (!tmp.empty()) {
        uf.unite(i, j);
      }
    }
  }
  std::map<unsigned int, std::vector<unsigned int>> members;
  for (unsigned int i = 0; i < rings.size(); ++i) {
    members[uf.find(i)].push_back(i);
  }
  groups.clear();
  groups.reserve(members.size());
  for (const auto &kv : members) {
    groups.push_back(kv.second);
  }
  std::vector<AtomClassification> classification(mol.getNumAtoms());
  for (unsigned int gidx = 0; gidx < groups.size(); ++gidx) {
    for (const auto ringIdx : groups[gidx]) {
      for (const auto atomIdx : rings[ringIdx]) {
        classification[atomIdx].inRing = true;
        classification[atomIdx].ringGroup = gidx;
      }
    }
  }
  return classification;
}

std::vector<unsigned int> computeAtomToVertex(const FeatTreeGraph &graph,
                                              unsigned int atomCount) {
  std::vector<unsigned int> mapping(atomCount,
                                    std::numeric_limits<unsigned int>::max());
  const auto nodeMap = boost::get(FeatTreeNode_t(), graph);
  for (auto vp = boost::vertices(graph); vp.first != vp.second; ++vp.first) {
    const auto &node = nodeMap[*vp.first];
    for (const auto atomIdx : node.atoms) {
      if (atomIdx < atomCount) {
        mapping[atomIdx] = static_cast<unsigned int>(*vp.first);
      }
    }
  }
  return mapping;
}

}  // namespace

void addRingsAndConnectors(const ROMol &mol, FeatTreeGraph &featGraph,
                           const FeatTreeParams &params) {
  std::vector<std::vector<int>> rings;
  std::vector<std::vector<unsigned int>> ringGroups;
  std::vector<unsigned int> ringSizes;
  const auto classification =
      classifyAtoms(mol, rings, ringGroups, ringSizes);

  std::vector<unsigned int> ringGroupToVertex(ringGroups.size(),
                                              std::numeric_limits<unsigned int>::max());
  auto nodeMap = boost::get(FeatTreeNode_t(), featGraph);

  for (unsigned int groupIdx = 0; groupIdx < ringGroups.size(); ++groupIdx) {
    std::vector<unsigned int> atoms;
    unsigned int minRing = 0u;
    unsigned int maxRing = 0u;
    for (const auto ringIdx : ringGroups[groupIdx]) {
      if (minRing == 0u || ringSizes[ringIdx] < minRing) {
        minRing = ringSizes[ringIdx];
      }
      if (ringSizes[ringIdx] > maxRing) {
        maxRing = ringSizes[ringIdx];
      }
      atoms.insert(atoms.end(), rings[ringIdx].begin(), rings[ringIdx].end());
    }
    finalizeAtomVector(atoms);
    FeatTreeNodeData node;
    node.atoms = atoms;
    node.kind = ringGroups[groupIdx].size() > 1 ? FeatTreeNodeKind::FusedRingSystem
                                                : FeatTreeNodeKind::RingSystem;
    node.minRingSize = static_cast<unsigned char>(std::min(minRing, 255u));
    node.maxRingSize = static_cast<unsigned char>(std::min(maxRing, 255u));
    accumulateNodeStats(mol, node.atoms, node);
    const auto v = boost::add_vertex(featGraph);
    nodeMap[v] = node;
    ringGroupToVertex[groupIdx] = static_cast<unsigned int>(v);
  }

  const unsigned int atomCount = mol.getNumAtoms();
  for (unsigned int idx = 0; idx < atomCount; ++idx) {
    const auto &ac = classification[idx];
    if (!ac.inRing) {
      const Atom *atom = mol.getAtomWithIdx(idx);
      PRECONDITION(atom, "missing atom during feature tree construction");
      const unsigned int deg = atom->getDegree() + atom->getTotalNumHs();
      if (deg > 1) {
        FeatTreeNodeData node;
        node.kind = FeatTreeNodeKind::Connector;
        node.atoms.push_back(idx);
        finalizeAtomVector(node.atoms);
        accumulateNodeStats(mol, node.atoms, node);
        auto v = boost::add_vertex(featGraph);
        nodeMap[v] = node;
      }
    }
  }
}

void addRingRingBonds(const ROMol &mol, FeatTreeGraph &featGraph) {
  const unsigned int atomCount = mol.getNumAtoms();
  auto atomToVertex = computeAtomToVertex(featGraph, atomCount);
  auto edgeMap = boost::get(FeatTreeEdge_t(), featGraph);
  for (const auto bond : mol.bonds()) {
    const auto beginIdx = bond->getBeginAtomIdx();
    const auto endIdx = bond->getEndAtomIdx();
    const auto v1 = atomToVertex[beginIdx];
    const auto v2 = atomToVertex[endIdx];
    if (v1 == std::numeric_limits<unsigned int>::max() ||
        v2 == std::numeric_limits<unsigned int>::max()) {
      continue;
    }
    if (v1 == v2) {
      continue;
    }
    FeatTreeEdgeData edgeData;
    const auto nodeMap = boost::get(FeatTreeNode_t(), featGraph);
    edgeData.ringEndCount = static_cast<unsigned char>(
        (isRingNode(nodeMap[v1]) ? 1 : 0) + (isRingNode(nodeMap[v2]) ? 1 : 0));
    auto res = boost::add_edge(v1, v2, featGraph);
    if (res.second) {
      edgeMap[res.first] = edgeData;
    }
  }
}

std::vector<unsigned int> addNonringAtoms(const ROMol &mol,
                                          FeatTreeGraph &featGraph,
                                          const FeatTreeParams &params) {
  (void)params;
  return computeAtomToVertex(featGraph, mol.getNumAtoms());
}

void addBondsFromNonringAtoms(const ROMol &mol, FeatTreeGraph &featGraph,
                              const std::vector<unsigned int> &atomIndices) {
  auto edgeMap = boost::get(FeatTreeEdge_t(), featGraph);
  const auto nodeMap = boost::get(FeatTreeNode_t(), featGraph);
  for (const auto bond : mol.bonds()) {
    const auto beginIdx = bond->getBeginAtomIdx();
    const auto endIdx = bond->getEndAtomIdx();
    if (beginIdx >= atomIndices.size() || endIdx >= atomIndices.size()) {
      continue;
    }
    const auto v1 = atomIndices[beginIdx];
    const auto v2 = atomIndices[endIdx];
    if (v1 == std::numeric_limits<unsigned int>::max() ||
        v2 == std::numeric_limits<unsigned int>::max() || v1 == v2) {
      continue;
    }
    FeatTreeEdgeData data;
    data.ringEndCount = static_cast<unsigned char>(
        (isRingNode(nodeMap[v1]) ? 1 : 0) + (isRingNode(nodeMap[v2]) ? 1 : 0));
    auto res = boost::add_edge(v1, v2, featGraph);
    if (res.second) {
      edgeMap[res.first] = data;
    }
  }
}

void addZeroNodes(FeatTreeGraph &featGraph, const FeatTreeParams &params) {
  if (!params.includeZeroNodes) {
    return;
  }
  auto nodeMap = boost::get(FeatTreeNode_t(), featGraph);
  auto edgeMap = boost::get(FeatTreeEdge_t(), featGraph);
  std::vector<std::pair<unsigned int, unsigned int>> pending;
  std::vector<boost::graph_traits<FeatTreeGraph>::edge_descriptor> toRemove;
  for (auto ep = boost::edges(featGraph); ep.first != ep.second; ++ep.first) {
    const auto src = boost::source(*ep.first, featGraph);
    const auto dst = boost::target(*ep.first, featGraph);
    if (!isRingNode(nodeMap[src]) && !isRingNode(nodeMap[dst])) {
      pending.emplace_back(src, dst);
      toRemove.push_back(*ep.first);
    }
  }
  for (const auto &edge : pending) {
    FeatTreeNodeData node;
    node.kind = FeatTreeNodeKind::ZeroNode;
    auto newV = boost::add_vertex(featGraph);
    nodeMap[newV] = node;
    FeatTreeEdgeData data;
    data.ringEndCount = 0u;
    auto e1 = boost::add_edge(edge.first, newV, featGraph);
    auto e2 = boost::add_edge(newV, edge.second, featGraph);
    if (e1.second) {
      edgeMap[e1.first] = data;
      setEdgeFlag(edgeMap[e1.first], FeatTreeEdgeFlags::PathCompressed, true);
    }
    if (e2.second) {
      edgeMap[e2.first] = data;
      setEdgeFlag(edgeMap[e2.first], FeatTreeEdgeFlags::PathCompressed, true);
    }
  }
  for (auto edge : toRemove) {
    boost::remove_edge(edge, featGraph);
  }
}

void replaceCycles(FeatTreeGraph &featGraph) {
  // fused ring cycles are already collapsed during construction; ensure
  // canonical numbering for determinism.
  canonicalizeFeatTree(featGraph);
}

}  // namespace FeatTrees
}  // namespace RDKit

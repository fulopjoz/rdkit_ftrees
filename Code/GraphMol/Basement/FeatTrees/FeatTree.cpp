//
//  Copyright (C) 2024  RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "FeatTree.h"

#include <GraphMol/RDKitBase.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/RingInfo.h>
#include <GraphMol/ChemicalFeatures/MolChemicalFeatureFactory.h>
#include <GraphMol/ChemicalFeatures/ChemicalFeature.h>
#include <RDGeneral/RDLog.h>
#include <RDConfig/Dirs.h>

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <unordered_map>

#include "FeatTreeUtils.h"

namespace RDKit {
namespace FeatTrees {
namespace {

using Vertex = boost::graph_traits<FeatTreeGraph>::vertex_descriptor;
using Edge = boost::graph_traits<FeatTreeGraph>::edge_descriptor;

unsigned int getDegree(const FeatTreeGraph &graph, Vertex v) {
  return static_cast<unsigned int>(boost::degree(v, graph));
}

const char *kindToString(FeatTreeNodeKind kind) {
  switch (kind) {
    case FeatTreeNodeKind::RingSystem:
      return "RingSystem";
    case FeatTreeNodeKind::FusedRingSystem:
      return "FusedRingSystem";
    case FeatTreeNodeKind::Connector:
      return "Connector";
    case FeatTreeNodeKind::BranchGroup:
      return "BranchGroup";
    case FeatTreeNodeKind::ZeroNode:
      return "ZeroNode";
    case FeatTreeNodeKind::FeatureGroup:
      return "FeatureGroup";
  }
  return "Unknown";
}

unsigned int kindPriority(FeatTreeNodeKind kind) {
  switch (kind) {
    case FeatTreeNodeKind::FusedRingSystem:
      return 0;
    case FeatTreeNodeKind::RingSystem:
      return 1;
    case FeatTreeNodeKind::Connector:
      return 2;
    case FeatTreeNodeKind::BranchGroup:
      return 3;
    case FeatTreeNodeKind::FeatureGroup:
      return 4;
    case FeatTreeNodeKind::ZeroNode:
      return 5;
  }
  return 6;
}

uint64_t hashNode(const FeatTreeNodeData &node) {
  uint64_t res = static_cast<uint64_t>(kindPriority(node.kind)) << 48;
  res |= static_cast<uint64_t>(node.atoms.empty() ? 0 : node.atoms.front()) << 32;
  res ^= static_cast<uint64_t>(node.aromaticAtomCount) << 24;
  res ^= static_cast<uint64_t>(node.heteroAtomCount) << 16;
  res ^= static_cast<uint64_t>(node.flags) << 8;
  res ^= static_cast<uint64_t>(node.maxRingSize);
  return res;
}

void mergeNodeData(FeatTreeNodeData &dst, const FeatTreeNodeData &src) {
  dst.atoms.insert(dst.atoms.end(), src.atoms.begin(), src.atoms.end());
  finalizeAtomVector(dst.atoms);
  dst.aromaticAtomCount =
      static_cast<unsigned char>(std::min<unsigned int>(
          static_cast<unsigned int>(dst.aromaticAtomCount) +
              static_cast<unsigned int>(src.aromaticAtomCount),
          255u));
  dst.heteroAtomCount =
      static_cast<unsigned char>(std::min<unsigned int>(
          static_cast<unsigned int>(dst.heteroAtomCount) +
              static_cast<unsigned int>(src.heteroAtomCount),
          255u));
  dst.flags |= src.flags;
  if (dst.minRingSize == 0u ||
      (src.minRingSize != 0u && src.minRingSize < dst.minRingSize)) {
    dst.minRingSize = src.minRingSize;
  }
  if (src.maxRingSize > dst.maxRingSize) {
    dst.maxRingSize = src.maxRingSize;
  }
}

std::vector<unsigned int> atomToVertex(const FeatTreeGraph &graph,
                                       unsigned int atomCount) {
  std::vector<unsigned int> mapping(atomCount,
                                    std::numeric_limits<unsigned int>::max());
  const auto nodeMap = boost::get(FeatTreeNode_t(), graph);
  for (auto vp = boost::vertices(graph); vp.first != vp.second; ++vp.first) {
    for (const auto atomIdx : nodeMap[*vp.first].atoms) {
      if (atomIdx < atomCount) {
        mapping[atomIdx] = static_cast<unsigned int>(*vp.first);
      }
    }
  }
  return mapping;
}

bool isConnector(const FeatTreeNodeData &node) {
  return node.kind == FeatTreeNodeKind::Connector;
}

void ensureInvariant(const FeatTreeGraph &graph) {
  const auto nodeMap = boost::get(FeatTreeNode_t(), graph);
  std::set<unsigned int> seenAtoms;
  for (auto vp = boost::vertices(graph); vp.first != vp.second; ++vp.first) {
    const auto &node = nodeMap[*vp.first];
    if (!isZeroNode(node)) {
      PRECONDITION(!node.atoms.empty(), "node with empty atom list");
    }
    for (auto idx : node.atoms) {
      PRECONDITION(seenAtoms.insert(idx).second ||
                       node.kind == FeatTreeNodeKind::FeatureGroup,
                   "atoms assigned to multiple structural nodes");
    }
  }
}

struct ChainResult {
  bool compressible = false;
  std::vector<Vertex> members;
  std::vector<Vertex> endpoints;
};

ChainResult gatherConnectorChain(const FeatTreeGraph &graph, Vertex start,
                                 const std::vector<bool> &eligible) {
  ChainResult result;
  result.members.push_back(start);
  auto neighbours = boost::adjacent_vertices(start, graph);
  if (std::distance(neighbours.first, neighbours.second) != 2) {
    result.compressible = false;
    return result;
  }
  result.compressible = true;
  std::set<Vertex> visited;
  visited.insert(start);
  for (Vertex nbr : {*(neighbours.first), *(std::next(neighbours.first))}) {
    Vertex prev = start;
    Vertex current = nbr;
    while (eligible[current]) {
      if (!visited.insert(current).second) {
        result.compressible = false;
        break;
      }
      result.members.push_back(current);
      auto adj = boost::adjacent_vertices(current, graph);
      Vertex next = *(adj.first) == prev ? *(std::next(adj.first)) : *(adj.first);
      prev = current;
      current = next;
    }
    result.endpoints.push_back(current);
  }
  if (result.endpoints.size() == 2 && result.endpoints[0] != result.endpoints[1] &&
      result.compressible) {
    return result;
  }
  result.compressible = false;
  return result;
}

void compressPaths(FeatTreeGraph &graph, const FeatTreeParams &params) {
  if (!params.compressPaths) {
    return;
  }
  auto nodeMap = boost::get(FeatTreeNode_t(), graph);
  std::vector<Vertex> vertices;
  vertices.reserve(boost::num_vertices(graph));
  for (auto vp = boost::vertices(graph); vp.first != vp.second; ++vp.first) {
    vertices.push_back(*vp.first);
  }
  std::vector<bool> eligible(boost::num_vertices(graph), false);
  for (auto v : vertices) {
    const auto &node = nodeMap[v];
    if (isConnector(node) && getDegree(graph, v) == 2u) {
      eligible[v] = true;
    }
  }
  std::vector<unsigned int> mapping(boost::num_vertices(graph),
                                    std::numeric_limits<unsigned int>::max());
  std::vector<bool> consumed(boost::num_vertices(graph), false);
  std::vector<bool> wasCompressed(boost::num_vertices(graph), false);
  std::vector<FeatTreeNodeData> newNodes;
  newNodes.reserve(vertices.size());

  for (auto v : vertices) {
    if (consumed[v]) {
      continue;
    }
    if (eligible[v]) {
      auto chain = gatherConnectorChain(graph, v, eligible);
      if (chain.compressible) {
        FeatTreeNodeData merged = nodeMap[chain.members.front()];
        merged.kind = FeatTreeNodeKind::Connector;
        for (std::size_t i = 1; i < chain.members.size(); ++i) {
          mergeNodeData(merged, nodeMap[chain.members[i]]);
        }
        const unsigned int newIdx = static_cast<unsigned int>(newNodes.size());
        newNodes.push_back(merged);
        for (auto member : chain.members) {
          consumed[member] = true;
          mapping[member] = newIdx;
          wasCompressed[member] = true;
        }
        continue;
      }
    }
    const unsigned int newIdx = static_cast<unsigned int>(newNodes.size());
    newNodes.push_back(nodeMap[v]);
    consumed[v] = true;
    mapping[v] = newIdx;
  }

  FeatTreeGraph rebuilt;
  auto rebuiltNodeMap = boost::get(FeatTreeNode_t(), rebuilt);
  for (const auto &node : newNodes) {
    auto nv = boost::add_vertex(rebuilt);
    rebuiltNodeMap[nv] = node;
  }
  std::map<std::pair<unsigned int, unsigned int>, FeatTreeEdgeData> edgeCache;
  auto edgeMap = boost::get(FeatTreeEdge_t(), graph);
  for (auto ep = boost::edges(graph); ep.first != ep.second; ++ep.first) {
    const auto src = boost::source(*ep.first, graph);
    const auto dst = boost::target(*ep.first, graph);
    const auto newSrc = mapping[src];
    const auto newDst = mapping[dst];
    if (newSrc == std::numeric_limits<unsigned int>::max() ||
        newDst == std::numeric_limits<unsigned int>::max() || newSrc == newDst) {
      continue;
    }
    auto key = std::minmax(newSrc, newDst);
    auto data = edgeMap[*ep.first];
    if (wasCompressed[src] || wasCompressed[dst]) {
      data.flags |= FeatTreeEdgeFlags::PathCompressed;
    }
    auto it = edgeCache.find(key);
    if (it == edgeCache.end()) {
      edgeCache[key] = data;
    } else {
      it->second.flags |= data.flags;
      it->second.ringEndCount = std::max(it->second.ringEndCount, data.ringEndCount);
    }
  }
  auto rebuiltEdgeMap = boost::get(FeatTreeEdge_t(), rebuilt);
  for (const auto &kv : edgeCache) {
    auto res = boost::add_edge(kv.first.first, kv.first.second, rebuilt);
    rebuiltEdgeMap[res.first] = kv.second;
  }
  graph = std::move(rebuilt);
}

void mergeBranchGroups(FeatTreeGraph &graph, const FeatTreeParams &params) {
  if (!params.mergeBranchGroups) {
    return;
  }
  auto nodeMap = boost::get(FeatTreeNode_t(), graph);
  std::vector<Vertex> vertices;
  for (auto vp = boost::vertices(graph); vp.first != vp.second; ++vp.first) {
    vertices.push_back(*vp.first);
  }
  for (auto v : vertices) {
    auto &node = nodeMap[v];
    if (node.kind == FeatTreeNodeKind::Connector &&
        getDegree(graph, v) == 1u &&
        node.atoms.size() <= params.maxBranchGroupSize) {
      node.kind = FeatTreeNodeKind::BranchGroup;
    }
  }
}

void annotateFeatureGroups(FeatTreeGraph &graph, const ROMol *mol,
                           const FeatTreeParams &params) {
  if (!params.annotateFeatures || !mol) {
    return;
  }
  std::unique_ptr<MolChemicalFeatureFactory> factory;
  try {
    std::string fdef = RDConfig::getDataDir() + "/BaseFeatures.fdef";
    factory.reset(new MolChemicalFeatureFactory(fdef));
  } catch (const std::exception &e) {
    BOOST_LOG(rdWarningLog)
        << "FeatTrees: unable to load BaseFeatures.fdef (" << e.what()
        << "), using heuristic feature annotation." << std::endl;
  }
  auto nodeMap = boost::get(FeatTreeNode_t(), graph);
  if (factory) {
    auto feats = factory->getFeaturesForMol(*mol);
    for (const auto &feat : feats) {
      FeatTreeNodeData node;
      node.kind = FeatTreeNodeKind::FeatureGroup;
      node.flags = 0u;
      node.atoms.assign(feat->getAtomIds().begin(), feat->getAtomIds().end());
      finalizeAtomVector(node.atoms);
      node.aromaticAtomCount = 0u;
      node.heteroAtomCount = static_cast<unsigned char>(node.atoms.size());
      if (feat->getFamily() == "Donor") {
        setNodeFlag(node, FeatTreeNodeFlags::Donor, true);
      }
      if (feat->getFamily() == "Acceptor") {
        setNodeFlag(node, FeatTreeNodeFlags::Acceptor, true);
      }
      auto v = boost::add_vertex(graph);
      nodeMap[v] = node;
    }
  } else {
    for (auto atom : mol->atoms()) {
      if (atom->getAtomicNum() == 0) {
        continue;
      }
      if (atom->getAtomicNum() != 6) {
        FeatTreeNodeData node;
        node.kind = FeatTreeNodeKind::FeatureGroup;
        node.atoms.push_back(atom->getIdx());
        finalizeAtomVector(node.atoms);
        node.heteroAtomCount = 1u;
        node.aromaticAtomCount = atom->getIsAromatic() ? 1u : 0u;
        setNodeFlag(node, FeatTreeNodeFlags::HasHetero, true);
        setNodeFlag(node, FeatTreeNodeFlags::Donor,
                    atom->getTotalNumHs(true) > 0);
        setNodeFlag(node, FeatTreeNodeFlags::Acceptor,
                    atom->getTotalValence() - atom->getExplicitValence() > 0);
        auto v = boost::add_vertex(graph);
        nodeMap[v] = node;
      }
    }
  }
}

struct NodeSignature {
  FeatTreeNodeKind kind;
  unsigned char atomCount;
  unsigned char aromaticCount;
  unsigned char heteroBin;
  unsigned char minRing;
  unsigned char maxRing;
  unsigned char flags;

  bool operator<(const NodeSignature &other) const {
    return std::tie(kind, atomCount, aromaticCount, heteroBin, minRing, maxRing,
                    flags) <
           std::tie(other.kind, other.atomCount, other.aromaticCount,
                    other.heteroBin, other.minRing, other.maxRing, other.flags);
  }
};

NodeSignature makeSignature(const FeatTreeNodeData &node) {
  NodeSignature sig;
  sig.kind = node.kind;
  sig.atomCount = static_cast<unsigned char>(
      std::min<std::size_t>(node.atoms.size(), static_cast<std::size_t>(255)));
  sig.aromaticCount = node.aromaticAtomCount;
  sig.heteroBin = static_cast<unsigned char>(
      std::min<unsigned int>(node.heteroAtomCount, 4u));
  sig.minRing = node.minRingSize;
  sig.maxRing = node.maxRingSize;
  sig.flags = node.flags;
  return sig;
}

double nodeWeight(const FeatTreeNodeData &node, const FeatTreeParams &params) {
  double base = 1.0 + std::log1p(static_cast<double>(node.atoms.size()));
  switch (node.kind) {
    case FeatTreeNodeKind::RingSystem:
    case FeatTreeNodeKind::FusedRingSystem:
      return params.ringWeight * base;
    case FeatTreeNodeKind::Connector:
    case FeatTreeNodeKind::BranchGroup:
      return params.connectorWeight * base;
    case FeatTreeNodeKind::FeatureGroup:
      return params.featureGroupWeight * base;
    case FeatTreeNodeKind::ZeroNode:
      return 0.5 * base;
  }
  return base;
}

std::map<NodeSignature, double> signatureWeights(const FeatTreeGraph &graph,
                                                 const FeatTreeParams &params) {
  std::map<NodeSignature, double> weights;
  const auto nodeMap = boost::get(FeatTreeNode_t(), graph);
  for (auto vp = boost::vertices(graph); vp.first != vp.second; ++vp.first) {
    const auto &node = nodeMap[*vp.first];
    if (node.kind == FeatTreeNodeKind::ZeroNode) {
      continue;
    }
    auto sig = makeSignature(node);
    weights[sig] += nodeWeight(node, params);
  }
  return weights;
}

}  // namespace

FeatTreeGraphSPtr molToBaseTree(const ROMol &mol, const FeatTreeParams &params) {
  auto result = FeatTreeGraphSPtr(new FeatTreeGraph());
  addRingsAndConnectors(mol, *result, params);
  addRingRingBonds(mol, *result);
  auto atomIndices = addNonringAtoms(mol, *result, params);
  addBondsFromNonringAtoms(mol, *result, atomIndices);
  addZeroNodes(*result, params);
  if (params.canonicalize) {
    canonicalizeFeatTree(*result);
  }
  return result;
}

void baseTreeToFeatTree(FeatTreeGraph &baseTree, const FeatTreeParams &params,
                        const ROMol *mol) {
  compressPaths(baseTree, params);
  mergeBranchGroups(baseTree, params);
  annotateFeatureGroups(baseTree, mol, params);
  if (params.canonicalize) {
    canonicalizeFeatTree(baseTree);
  }
  ensureInvariant(baseTree);
}

FeatTreeGraphSPtr molToFeatTree(const ROMol &mol, const FeatTreeParams &params) {
  auto tree = molToBaseTree(mol, params);
  baseTreeToFeatTree(*tree, params, &mol);
  return tree;
}

void canonicalizeFeatTree(FeatTreeGraph &graph) {
  auto nodeMap = boost::get(FeatTreeNode_t(), graph);
  std::vector<Vertex> vertices;
  for (auto vp = boost::vertices(graph); vp.first != vp.second; ++vp.first) {
    vertices.push_back(*vp.first);
  }
  struct CanonKey {
    Vertex v;
    uint64_t hash;
    unsigned int priority;
    unsigned int smallestAtom;
    std::size_t size;
  };
  std::vector<CanonKey> keys;
  keys.reserve(vertices.size());
  for (auto v : vertices) {
    const auto &node = nodeMap[v];
    CanonKey key;
    key.v = v;
    key.hash = hashNode(node);
    key.priority = kindPriority(node.kind);
    key.size = node.atoms.size();
    key.smallestAtom = node.atoms.empty() ? std::numeric_limits<unsigned int>::max()
                                          : node.atoms.front();
    keys.push_back(key);
  }
  std::sort(keys.begin(), keys.end(), [](const CanonKey &lhs, const CanonKey &rhs) {
    return std::tie(lhs.priority, lhs.smallestAtom, lhs.size, lhs.hash) <
           std::tie(rhs.priority, rhs.smallestAtom, rhs.size, rhs.hash);
  });
  std::vector<unsigned int> mapping(boost::num_vertices(graph));
  FeatTreeGraph reordered;
  auto reorderedNodeMap = boost::get(FeatTreeNode_t(), reordered);
  for (unsigned int newIdx = 0; newIdx < keys.size(); ++newIdx) {
    mapping[keys[newIdx].v] = newIdx;
    auto nv = boost::add_vertex(reordered);
    reorderedNodeMap[nv] = nodeMap[keys[newIdx].v];
  }
  auto edgeMap = boost::get(FeatTreeEdge_t(), graph);
  std::map<std::pair<unsigned int, unsigned int>, FeatTreeEdgeData> edgeCache;
  for (auto ep = boost::edges(graph); ep.first != ep.second; ++ep.first) {
    auto src = mapping[boost::source(*ep.first, graph)];
    auto dst = mapping[boost::target(*ep.first, graph)];
    if (src == dst) {
      continue;
    }
    auto key = std::minmax(src, dst);
    auto data = edgeMap[*ep.first];
    auto it = edgeCache.find(key);
    if (it == edgeCache.end()) {
      edgeCache[key] = data;
    } else {
      it->second.flags |= data.flags;
      it->second.ringEndCount = std::max(it->second.ringEndCount, data.ringEndCount);
    }
  }
  auto reorderedEdgeMap = boost::get(FeatTreeEdge_t(), reordered);
  for (const auto &kv : edgeCache) {
    auto res = boost::add_edge(kv.first.first, kv.first.second, reordered);
    reorderedEdgeMap[res.first] = kv.second;
  }
  graph = std::move(reordered);
}

std::string featTreeToJSON(const FeatTreeGraph &graph) {
  std::ostringstream oss;
  oss << "{\"nodes\":[";
  const auto nodeMap = boost::get(FeatTreeNode_t(), graph);
  bool first = true;
  unsigned int idx = 0;
  for (auto vp = boost::vertices(graph); vp.first != vp.second; ++vp.first, ++idx) {
    if (!first) {
      oss << ',';
    }
    first = false;
    const auto &node = nodeMap[*vp.first];
    oss << "{\"index\":" << idx;
    oss << ",\"kind\":\"" << kindToString(node.kind) << "\"";
    oss << ",\"atoms\":[";
    for (std::size_t i = 0; i < node.atoms.size(); ++i) {
      if (i) {
        oss << ',';
      }
      oss << node.atoms[i];
    }
    oss << "]";
    oss << ",\"flags\":" << static_cast<unsigned int>(node.flags);
    oss << ",\"aromaticAtomCount\":" << static_cast<unsigned int>(node.aromaticAtomCount);
    oss << ",\"heteroAtomCount\":" << static_cast<unsigned int>(node.heteroAtomCount);
    oss << ",\"minRingSize\":" << static_cast<unsigned int>(node.minRingSize);
    oss << ",\"maxRingSize\":" << static_cast<unsigned int>(node.maxRingSize);
    oss << '}';
  }
  oss << "],\"edges\":[";
  const auto edgeMap = boost::get(FeatTreeEdge_t(), graph);
  first = true;
  for (auto ep = boost::edges(graph); ep.first != ep.second; ++ep.first) {
    if (!first) {
      oss << ',';
    }
    first = false;
    oss << "{\"source\":" << boost::source(*ep.first, graph);
    oss << ",\"target\":" << boost::target(*ep.first, graph);
    const auto &edge = edgeMap[*ep.first];
    oss << ",\"ringEndCount\":" << static_cast<unsigned int>(edge.ringEndCount);
    oss << ",\"flags\":" << static_cast<unsigned int>(edge.flags);
    oss << '}';
  }
  oss << "]}";
  return oss.str();
}

double calcFeatTreeSimilarity(const FeatTreeGraph &g1, const FeatTreeGraph &g2,
                              const FeatTreeParams &params) {
  auto weights1 = signatureWeights(g1, params);
  auto weights2 = signatureWeights(g2, params);
  auto it1 = weights1.begin();
  auto it2 = weights2.begin();
  double inter = 0.0;
  double uni = 0.0;
  while (it1 != weights1.end() && it2 != weights2.end()) {
    if (it1->first < it2->first) {
      uni += it1->second;
      ++it1;
    } else if (it2->first < it1->first) {
      uni += it2->second;
      ++it2;
    } else {
      inter += std::min(it1->second, it2->second);
      uni += std::max(it1->second, it2->second);
      ++it1;
      ++it2;
    }
  }
  for (; it1 != weights1.end(); ++it1) {
    uni += it1->second;
  }
  for (; it2 != weights2.end(); ++it2) {
    uni += it2->second;
  }
  if (uni <= 0.0) {
    return 0.0;
  }
  return inter / uni;
}

double calcFeatTreeSimilarity(const ROMol &mol1, const ROMol &mol2,
                              const FeatTreeParams &params) {
  auto tree1 = molToFeatTree(mol1, params);
  auto tree2 = molToFeatTree(mol2, params);
  return calcFeatTreeSimilarity(*tree1, *tree2, params);
}

}  // namespace FeatTrees
}  // namespace RDKit

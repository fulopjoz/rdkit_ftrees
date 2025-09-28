//
//  Copyright (C) 2024  RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/test.h>
#include <GraphMol/Basement/FeatTrees/FeatTree.h>
#include <GraphMol/Basement/FeatTrees/FeatTreeUtils.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <algorithm>
#include <memory>
#include <set>
#include <string>

using namespace RDKit;
using namespace RDKit::FeatTrees;

namespace {

std::set<std::pair<unsigned int, unsigned int>> edgePairs(const FeatTreeGraph &g) {
  std::set<std::pair<unsigned int, unsigned int>> res;
  for (auto ep = boost::edges(g); ep.first != ep.second; ++ep.first) {
    auto src = static_cast<unsigned int>(boost::source(*ep.first, g));
    auto dst = static_cast<unsigned int>(boost::target(*ep.first, g));
    if (src > dst) {
      std::swap(src, dst);
    }
    res.emplace(src, dst);
  }
  return res;
}

}  // namespace

void testSimpleChainBaseTree() {
  auto mol = std::unique_ptr<ROMol>(SmilesToMol("CCC"));
  TEST_ASSERT(mol);
  FeatTreeParams params;
  params.includeZeroNodes = false;
  auto baseTree = molToBaseTree(*mol, params);
  TEST_ASSERT(baseTree);
  TEST_ASSERT(boost::num_vertices(*baseTree) == 3);
  TEST_ASSERT(boost::num_edges(*baseTree) == 2);
  auto nodeMap = boost::get(FeatTreeNode_t(), *baseTree);
  for (auto vp = boost::vertices(*baseTree); vp.first != vp.second; ++vp.first) {
    const auto &node = nodeMap[*vp.first];
    TEST_ASSERT(node.kind == FeatTreeNodeKind::Connector);
    TEST_ASSERT(node.atoms.size() == 1);
  }
  auto json1 = featTreeToJSON(*baseTree);
  auto json2 = featTreeToJSON(*baseTree);
  TEST_ASSERT(json1 == json2);
}

void testFusedRingDetection() {
  auto mol = std::unique_ptr<ROMol>(SmilesToMol("C1CCC2CCCCC12"));
  TEST_ASSERT(mol);
  auto baseTree = molToBaseTree(*mol);
  TEST_ASSERT(baseTree);
  auto nodeMap = boost::get(FeatTreeNode_t(), *baseTree);
  unsigned int fusedCount = 0;
  for (auto vp = boost::vertices(*baseTree); vp.first != vp.second; ++vp.first) {
    if (nodeMap[*vp.first].kind == FeatTreeNodeKind::FusedRingSystem) {
      ++fusedCount;
    }
  }
  TEST_ASSERT(fusedCount == 1);
}

void testPathCompression() {
  auto mol = std::unique_ptr<ROMol>(SmilesToMol("CCCCCC"));
  TEST_ASSERT(mol);
  FeatTreeParams params;
  params.compressPaths = true;
  params.mergeBranchGroups = false;
  auto baseTree = molToBaseTree(*mol, params);
  TEST_ASSERT(boost::num_vertices(*baseTree) == 6);
  baseTreeToFeatTree(*baseTree, params, mol.get());
  TEST_ASSERT(boost::num_vertices(*baseTree) <= 4);
  unsigned int branchCount = 0;
  auto nodeMap = boost::get(FeatTreeNode_t(), *baseTree);
  for (auto vp = boost::vertices(*baseTree); vp.first != vp.second; ++vp.first) {
    if (nodeMap[*vp.first].kind == FeatTreeNodeKind::Connector) {
      ++branchCount;
    }
  }
  TEST_ASSERT(branchCount <= 4);
}

void testBranchGrouping() {
  auto mol = std::unique_ptr<ROMol>(SmilesToMol("CC(C)C"));
  TEST_ASSERT(mol);
  FeatTreeParams params;
  params.maxBranchGroupSize = 1;
  params.mergeBranchGroups = true;
  auto tree = molToFeatTree(*mol, params);
  auto nodeMap = boost::get(FeatTreeNode_t(), *tree);
  bool seenBranchGroup = false;
  for (auto vp = boost::vertices(*tree); vp.first != vp.second; ++vp.first) {
    if (nodeMap[*vp.first].kind == FeatTreeNodeKind::BranchGroup) {
      seenBranchGroup = true;
    }
  }
  TEST_ASSERT(seenBranchGroup);
}

void testFeatureAnnotationFallback() {
  auto mol = std::unique_ptr<ROMol>(SmilesToMol("CCO"));
  TEST_ASSERT(mol);
  FeatTreeParams params;
  params.annotateFeatures = true;
  params.mergeBranchGroups = false;
  auto tree = molToFeatTree(*mol, params);
  auto nodeMap = boost::get(FeatTreeNode_t(), *tree);
  unsigned int featureCount = 0;
  for (auto vp = boost::vertices(*tree); vp.first != vp.second; ++vp.first) {
    if (nodeMap[*vp.first].kind == FeatTreeNodeKind::FeatureGroup) {
      ++featureCount;
    }
  }
  TEST_ASSERT(featureCount >= 1);
}

void testDeterminism() {
  auto mol = std::unique_ptr<ROMol>(SmilesToMol("C1=CC=CC=C1"));
  TEST_ASSERT(mol);
  FeatTreeParams params;
  params.canonicalize = true;
  auto tree1 = molToFeatTree(*mol, params);
  auto tree2 = molToFeatTree(*mol, params);
  TEST_ASSERT(featTreeToJSON(*tree1) == featTreeToJSON(*tree2));
  TEST_ASSERT(edgePairs(*tree1) == edgePairs(*tree2));
}

int main() {
  RDLog::InitLogs();
  testSimpleChainBaseTree();
  testFusedRingDetection();
  testPathCompression();
  testBranchGrouping();
  testFeatureAnnotationFallback();
  testDeterminism();
  return 0;
}

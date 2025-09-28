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
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <RDGeneral/Invariant.h>

#include <functional>
#include <memory>

using namespace RDKit;
using namespace RDKit::FeatTrees;

namespace {

bool throwsInvariant(const std::function<void()> &fn) {
  try {
    fn();
  } catch (const Invariant::InvariantViolationException &) {
    return true;
  }
  return false;
}

}  // namespace

void testValidateParamsFailures() {
  FeatTreeParams params;
  params.ringWeight = -1.0;
  TEST_ASSERT(throwsInvariant([&]() { validateParams(params); }));
  params.ringWeight = 1.0;
  params.connectorWeight = 0.0;
  TEST_ASSERT(throwsInvariant([&]() { validateParams(params); }));
  params.connectorWeight = 1.0;
  params.featureGroupWeight = -0.1;
  TEST_ASSERT(throwsInvariant([&]() { validateParams(params); }));
  params.featureGroupWeight = 1.0;
  params.maxBranchGroupSize = 20;
  TEST_ASSERT(throwsInvariant([&]() { validateParams(params); }));
  params.maxBranchGroupSize = 4;
  params.similarityAutoThreshold = 0;
  TEST_ASSERT(throwsInvariant([&]() { validateParams(params); }));
}

void testDuplicateAtomDetection() {
  FeatTreeGraph graph;
  auto nodeMap = boost::get(FeatTreeNode_t(), graph);
  auto v = boost::add_vertex(graph);
  nodeMap[v].kind = FeatTreeNodeKind::Connector;
  nodeMap[v].atoms = {0, 0};
  FeatTreeParams params;
  TEST_ASSERT(throwsInvariant([&]() { validateFeatTree(graph, params, true); }));
}

void testZeroNodePolicy() {
  FeatTreeGraph graph;
  auto nodeMap = boost::get(FeatTreeNode_t(), graph);
  auto v = boost::add_vertex(graph);
  nodeMap[v].kind = FeatTreeNodeKind::ZeroNode;
  FeatTreeParams params;
  params.includeZeroNodes = false;
  TEST_ASSERT(throwsInvariant([&]() { validateFeatTree(graph, params, true); }));
  params.includeZeroNodes = true;
  validateFeatTree(graph, params, true);
}

void testAllowZeroNodesFlag() {
  FeatTreeGraph graph;
  FeatTreeParams params;
  TEST_ASSERT(throwsInvariant([&]() { validateFeatTree(graph, params, false); }));
  validateFeatTree(graph, params, true);
}

void testCanonicalOrderCheck() {
  auto mol = std::unique_ptr<ROMol>(SmilesToMol("c1ccccc1"));
  TEST_ASSERT(mol);
  FeatTreeParams params;
  params.canonicalize = true;
  auto tree = molToFeatTree(*mol, params);
  auto nodeMap = boost::get(FeatTreeNode_t(), *tree);
  if (boost::num_vertices(*tree) > 1) {
    std::swap(nodeMap[0], nodeMap[1]);
    TEST_ASSERT(throwsInvariant(
        [&]() { validateFeatTree(*tree, params, /*allowZeroNodes=*/true); }));
  }
}

void testMoleculeCoverage() {
  FeatTreeParams params;
  params.canonicalize = true;

  ROMol emptyMol;
  auto emptyTree = molToFeatTree(emptyMol, params);
  validateFeatTree(*emptyTree, params, true);
  TEST_ASSERT(boost::num_vertices(*emptyTree) == 0);

  auto macro = std::unique_ptr<ROMol>(SmilesToMol("C1CCCCCCCCCCCC1"));
  TEST_ASSERT(macro);
  auto macroTree = molToFeatTree(*macro, params);
  validateFeatTree(*macroTree, params, true);

  auto indole = std::unique_ptr<ROMol>(SmilesToMol("c1ccc2[nH]ccc2c1"));
  TEST_ASSERT(indole);
  auto indoleTree = molToFeatTree(*indole, params);
  validateFeatTree(*indoleTree, params, true);

  auto quinazoline =
      std::unique_ptr<ROMol>(SmilesToMol("c1nc2ncccc2n1"));
  TEST_ASSERT(quinazoline);
  auto quinTree = molToFeatTree(*quinazoline, params);
  validateFeatTree(*quinTree, params, true);

  auto spiro = std::unique_ptr<ROMol>(SmilesToMol("C1CCC2(CC1)CCC2"));
  TEST_ASSERT(spiro);
  auto spiroTree = molToFeatTree(*spiro, params);
  validateFeatTree(*spiroTree, params, true);

  auto salt = std::unique_ptr<ROMol>(SmilesToMol("CC.CN"));
  TEST_ASSERT(salt);
  auto saltTree = molToFeatTree(*salt, params);
  validateFeatTree(*saltTree, params, true);
  TEST_ASSERT(boost::num_vertices(*saltTree) >= 2);
}

void testBranchGroupToggleCount() {
  auto mol = std::unique_ptr<ROMol>(SmilesToMol("CC(C)C"));
  TEST_ASSERT(mol);
  FeatTreeParams params;
  params.mergeBranchGroups = false;
  auto withoutMerge = molToFeatTree(*mol, params);
  params.mergeBranchGroups = true;
  auto withMerge = molToFeatTree(*mol, params);
  TEST_ASSERT(boost::num_vertices(*withMerge) <=
              boost::num_vertices(*withoutMerge));
}

int main() {
  RDLog::InitLogs();
  testValidateParamsFailures();
  testDuplicateAtomDetection();
  testZeroNodePolicy();
  testAllowZeroNodesFlag();
  testCanonicalOrderCheck();
  testMoleculeCoverage();
  testBranchGroupToggleCount();
  return 0;
}

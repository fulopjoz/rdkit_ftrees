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
#include <GraphMol/Basement/FeatTrees/FeatTreeSimilarity.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <cmath>
#include <memory>

using namespace RDKit;
using namespace RDKit::FeatTrees;

void testEmptyGraphSimilarity() {
  FeatTreeGraph g1;
  FeatTreeGraph g2;
  FeatTreeGraph g3;
  auto nodeMap = boost::get(FeatTreeNode_t(), g3);
  auto v = boost::add_vertex(g3);
  nodeMap[v].kind = FeatTreeNodeKind::Connector;
  nodeMap[v].atoms = {0};
  FeatTreeParams params;
  const double simEmpty =
      calcFeatTreeSimilarity(g1, g2, FeatTreeSimilarityMethod::WeightedJaccard,
                             params);
  TEST_ASSERT(std::abs(simEmpty - 1.0) < 1e-8);
  const double simMixed =
      calcFeatTreeSimilarity(g1, g3, FeatTreeSimilarityMethod::WeightedJaccard,
                             params);
  TEST_ASSERT(std::abs(simMixed) < 1e-8);
}

void testIdenticalMoleculesSimilarity() {
  auto mol = std::unique_ptr<ROMol>(SmilesToMol("c1ccccc1"));
  TEST_ASSERT(mol);
  FeatTreeParams params;
  const double simWeighted = calcFeatTreeSimilarity(
      *mol, *mol, FeatTreeSimilarityMethod::WeightedJaccard, params);
  TEST_ASSERT(std::abs(simWeighted - 1.0) < 1e-8);
  const double simApprox = calcFeatTreeSimilarity(
      *mol, *mol, FeatTreeSimilarityMethod::ApproxEdit, params);
  TEST_ASSERT(std::abs(simApprox - 1.0) < 1e-8);
}

void testSubstructureMonotonicity() {
  auto benzene = std::unique_ptr<ROMol>(SmilesToMol("c1ccccc1"));
  auto toluene = std::unique_ptr<ROMol>(SmilesToMol("Cc1ccccc1"));
  auto butane = std::unique_ptr<ROMol>(SmilesToMol("CCCC"));
  TEST_ASSERT(benzene && toluene && butane);
  FeatTreeParams params;
  const double simBT = calcFeatTreeSimilarity(*benzene, *toluene,
                                              FeatTreeSimilarityMethod::Auto,
                                              params);
  const double simBB = calcFeatTreeSimilarity(*benzene, *butane,
                                              FeatTreeSimilarityMethod::Auto,
                                              params);
  const double simDisjoint = calcFeatTreeSimilarity(*toluene, *butane,
                                                    FeatTreeSimilarityMethod::Auto,
                                                    params);
  TEST_ASSERT(simBT >= simBB);
  TEST_ASSERT(simBB >= simDisjoint);
}

void testAutoSelectsApproximateForSmallGraphs() {
  auto ethane = std::unique_ptr<ROMol>(SmilesToMol("CC"));
  auto ethanol = std::unique_ptr<ROMol>(SmilesToMol("CCO"));
  TEST_ASSERT(ethane && ethanol);
  FeatTreeParams params;
  params.similarityAutoThreshold = 64;
  const double simApprox = calcFeatTreeSimilarity(*ethane, *ethanol,
                                                  FeatTreeSimilarityMethod::ApproxEdit,
                                                  params);
  const double simAuto = calcFeatTreeSimilarity(*ethane, *ethanol,
                                                FeatTreeSimilarityMethod::Auto,
                                                params);
  TEST_ASSERT(std::abs(simApprox - simAuto) < 1e-6);
}

int main() {
  RDLog::InitLogs();
  testEmptyGraphSimilarity();
  testIdenticalMoleculesSimilarity();
  testSubstructureMonotonicity();
  testAutoSelectsApproximateForSmallGraphs();
  return 0;
}

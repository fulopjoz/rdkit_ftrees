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

#include <memory>

using namespace RDKit;
using namespace RDKit::FeatTrees;

void testSimilarityIdentity() {
  auto mol = std::unique_ptr<ROMol>(SmilesToMol("CCO"));
  TEST_ASSERT(mol);
  FeatTreeParams params;
  params.annotateFeatures = true;
  const double sim = calcFeatTreeSimilarity(*mol, *mol, params);
  TEST_ASSERT(feq(sim, 1.0));
}

void testSimilarityDissimilar() {
  auto mol1 = std::unique_ptr<ROMol>(SmilesToMol("c1ccccc1"));
  auto mol2 = std::unique_ptr<ROMol>(SmilesToMol("CCCCC"));
  TEST_ASSERT(mol1 && mol2);
  FeatTreeParams params;
  const double sim = calcFeatTreeSimilarity(*mol1, *mol2, params);
  TEST_ASSERT(sim < 1.0);
  TEST_ASSERT(sim >= 0.0);
}

void testEditDistanceApprox() {
  auto mol1 = std::unique_ptr<ROMol>(SmilesToMol("c1ccccc1"));
  auto mol2 = std::unique_ptr<ROMol>(SmilesToMol("c1ccncc1"));
  TEST_ASSERT(mol1 && mol2);
  FeatTreeParams params;
  const double dist = calcFeatTreeEditDistanceApprox(*mol1, *mol2, params);
  TEST_ASSERT(dist >= 0.0);
}

int main() {
  RDLog::InitLogs();
  testSimilarityIdentity();
  testSimilarityDissimilar();
  testEditDistanceApprox();
  return 0;
}

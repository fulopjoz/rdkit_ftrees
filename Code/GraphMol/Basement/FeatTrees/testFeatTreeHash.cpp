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

#include <memory>

using namespace RDKit;
using namespace RDKit::FeatTrees;

void testHashConsistencyAcrossBuilders() {
  auto mol = std::unique_ptr<ROMol>(SmilesToMol("c1ccccc1"));
  TEST_ASSERT(mol);
  FeatTreeParams params;
  params.canonicalize = true;
  auto tree = molToFeatTree(*mol, params);
  auto base = molToBaseTree(*mol, params);
  baseTreeToFeatTree(*base, params, mol.get());
  const auto hash1 = hashFeatTree(*tree);
  const auto hash2 = hashFeatTree(*base);
  TEST_ASSERT(hash1 == hash2);
}

void testCanonicalizationIdempotent() {
  auto mol = std::unique_ptr<ROMol>(SmilesToMol("CC(C)C"));
  TEST_ASSERT(mol);
  FeatTreeParams params;
  params.canonicalize = true;
  auto tree = molToFeatTree(*mol, params);
  const auto initialHash = hashFeatTree(*tree);
  canonicalizeFeatTree(*tree);
  const auto hashAfterCanon = hashFeatTree(*tree);
  canonicalizeFeatTree(*tree);
  const auto hashAfterSecond = hashFeatTree(*tree);
  TEST_ASSERT(initialHash == hashAfterCanon);
  TEST_ASSERT(hashAfterCanon == hashAfterSecond);
}

void testHashIgnoresVertexOrdering() {
  auto mol = std::unique_ptr<ROMol>(SmilesToMol("c1ccc(C)cc1"));
  TEST_ASSERT(mol);
  FeatTreeParams params;
  params.canonicalize = true;
  auto tree = molToFeatTree(*mol, params);
  auto permuted = FeatTreeGraph(*tree);
  if (boost::num_vertices(permuted) > 1) {
    auto nodeMap = boost::get(FeatTreeNode_t(), permuted);
    std::swap(nodeMap[0], nodeMap[boost::num_vertices(permuted) - 1]);
  }
  const auto hash1 = hashFeatTree(*tree);
  const auto hash2 = hashFeatTree(permuted);
  TEST_ASSERT(hash1 == hash2);
}

int main() {
  RDLog::InitLogs();
  testHashConsistencyAcrossBuilders();
  testCanonicalizationIdempotent();
  testHashIgnoresVertexOrdering();
  return 0;
}

//
//  Copyright (C) 2024  RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "FeatTreeSimilarity.h"

#include <GraphMol/RDKitBase.h>

namespace RDKit {
namespace FeatTrees {

namespace {
double graphSize(const FeatTreeGraph &g) {
  return static_cast<double>(boost::num_vertices(g));
}
}  // namespace

double calcFeatTreeEditDistanceApprox(const FeatTreeGraph &g1,
                                      const FeatTreeGraph &g2,
                                      const FeatTreeParams &params) {
  const double sim = calcFeatTreeSimilarity(g1, g2, params);
  const double avgSize = 0.5 * (graphSize(g1) + graphSize(g2));
  return (1.0 - sim) * avgSize;
}

double calcFeatTreeEditDistanceApprox(const ROMol &mol1, const ROMol &mol2,
                                      const FeatTreeParams &params) {
  auto tree1 = molToFeatTree(mol1, params);
  auto tree2 = molToFeatTree(mol2, params);
  return calcFeatTreeEditDistanceApprox(*tree1, *tree2, params);
}

}  // namespace FeatTrees
}  // namespace RDKit

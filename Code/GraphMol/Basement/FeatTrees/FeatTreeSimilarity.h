//
//  Copyright (C) 2024  RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

/*! \file FeatTreeSimilarity.h

    \brief Helpers that expose edit-distance approximations for feature trees.
*/

#ifndef RD_FEATTREE_SIMILARITY_H
#define RD_FEATTREE_SIMILARITY_H

#include <RDGeneral/export.h>

#include "FeatTree.h"

namespace RDKit {
namespace FeatTrees {

//! \brief Approximates a tree edit distance using feature tree similarity.
RDKIT_GRAPHMOL_EXPORT double calcFeatTreeEditDistanceApprox(
    const FeatTreeGraph &g1, const FeatTreeGraph &g2,
    const FeatTreeParams &params = FeatTreeParams());

//! \brief Convenience overload operating directly on molecules.
RDKIT_GRAPHMOL_EXPORT double calcFeatTreeEditDistanceApprox(
    const ROMol &mol1, const ROMol &mol2,
    const FeatTreeParams &params = FeatTreeParams());

}  // namespace FeatTrees
}  // namespace RDKit

#endif  // RD_FEATTREE_SIMILARITY_H

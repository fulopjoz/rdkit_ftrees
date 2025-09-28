//
//  Copyright (C) 2024  RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

/*! \file FeatTreeUtils.h

  \brief Helper routines that expose intermediate steps of the feature tree
         construction pipeline.

  These are primarily used from unit tests to exercise individual stages of the
  builder.  End users should generally call molToFeatTree().
*/

#ifndef RD_FEATTREEUTILS_H
#define RD_FEATTREEUTILS_H

#include <RDGeneral/export.h>

#include "FeatTree.h"

#include <map>
#include <vector>

namespace RDKit {
class ROMol;

namespace FeatTrees {

//! \brief Returns indices of rings grouped by fused ring membership.
typedef std::map<unsigned int, std::vector<unsigned int>> RingGroupMap;

RDKIT_GRAPHMOL_EXPORT void addRingsAndConnectors(const ROMol &mol,
                                                 FeatTreeGraph &featGraph,
                                                 const FeatTreeParams &params);

RDKIT_GRAPHMOL_EXPORT void addRingRingBonds(const ROMol &mol,
                                            FeatTreeGraph &featGraph);

RDKIT_GRAPHMOL_EXPORT std::vector<unsigned int> addNonringAtoms(
    const ROMol &mol, FeatTreeGraph &featGraph,
    const FeatTreeParams &params);

RDKIT_GRAPHMOL_EXPORT void addBondsFromNonringAtoms(
    const ROMol &mol, FeatTreeGraph &featGraph,
    const std::vector<unsigned int> &atomIndices);

RDKIT_GRAPHMOL_EXPORT void addZeroNodes(FeatTreeGraph &featGraph,
                                        const FeatTreeParams &params);

RDKIT_GRAPHMOL_EXPORT void replaceCycles(FeatTreeGraph &featGraph);

}  // namespace FeatTrees
}  // namespace RDKit

#endif  // RD_FEATTREEUTILS_H

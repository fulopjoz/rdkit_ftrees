//
//  Copyright (C) 2024  RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RD_FEATTREE_SIMILARITY_H
#define RD_FEATTREE_SIMILARITY_H

#include <GraphMol/RDKitBase.h>
#include <GraphMol/Basement/FeatTrees/FeatTree.h>

#include <set>
#include <vector>

namespace RDKit {
namespace FeatTrees {

struct FragmentProfile {
  int volume = 0;
  bool ring = false;
  int donors = 0;
  int acceptors = 0;
  bool amide = false;
  bool aromatic = false;
  double hydrophobicity = 0.0;
};

FragmentProfile computeFragmentProfile(const ROMol &mol,
                                       const std::set<unsigned int> &indices);

std::vector<FragmentProfile> buildFragmentProfiles(const ROMol &mol,
                                                   const FeatTreeGraph &tree);

std::vector<FragmentProfile> buildFragmentProfiles(const ROMol &mol);

double compareFragmentProfiles(const FragmentProfile &lhs,
                               const FragmentProfile &rhs);

double compareProfileSets(const std::vector<FragmentProfile> &lhs,
                          const std::vector<FragmentProfile> &rhs);

double compareMolecules(const ROMol &lhs, const ROMol &rhs);

}  // namespace FeatTrees
}  // namespace RDKit

#endif  // RD_FEATTREE_SIMILARITY_H


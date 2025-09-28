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

#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/RingInfo.h>

#include <algorithm>
#include <cmath>
#include <memory>

namespace RDKit {
namespace FeatTrees {
namespace {

bool isDonor(const Atom *atom) {
  if (!atom) {
    return false;
  }
  const unsigned int atomicNum = atom->getAtomicNum();
  if ((atomicNum == 7 || atomicNum == 8 || atomicNum == 16) &&
      atom->getTotalNumHs() > 0) {
    return true;
  }
  return false;
}

bool isAcceptor(const Atom *atom) {
  if (!atom) {
    return false;
  }
  const unsigned int atomicNum = atom->getAtomicNum();
  if (atomicNum == 7 || atomicNum == 8 || atomicNum == 9 || atomicNum == 16) {
    if (atom->getFormalCharge() <= 0) {
      return true;
    }
  }
  return false;
}

bool containsAmideLike(const ROMol &mol,
                       const std::set<unsigned int> &indices) {
  for (const unsigned int idx : indices) {
    const Atom *atom = mol.getAtomWithIdx(idx);
    if (!atom || atom->getAtomicNum() != 6) {
      continue;
    }
    for (const auto bond : mol.atomBonds(atom)) {
      const Bond *b = bond;
      const unsigned int otherIdx = b->getOtherAtomIdx(idx);
      if (!indices.count(otherIdx)) {
        continue;
      }
      const Atom *nbr = mol.getAtomWithIdx(otherIdx);
      if (b->getBondType() == Bond::BondType::DOUBLE &&
          nbr->getAtomicNum() == 8) {
        for (const auto bond2 : mol.atomBonds(atom)) {
          const Bond *b2 = bond2;
          const unsigned int other = b2->getOtherAtomIdx(idx);
          if (!indices.count(other)) {
            continue;
          }
          const Atom *nbr2 = mol.getAtomWithIdx(other);
          if (nbr2->getAtomicNum() == 7 &&
              b2->getBondType() == Bond::BondType::SINGLE) {
            return true;
          }
        }
      }
    }
  }
  return false;
}

}  // namespace

FragmentProfile computeFragmentProfile(const ROMol &mol,
                                       const std::set<unsigned int> &indices) {
  FragmentProfile profile;
  int carbonCount = 0;
  for (const unsigned int idx : indices) {
    const Atom *atom = mol.getAtomWithIdx(idx);
    if (!atom) {
      continue;
    }
    ++profile.volume;
    if (atom->getAtomicNum() == 6) {
      ++carbonCount;
    }
    if (mol.getRingInfo()->numAtomRings(idx) > 0) {
      profile.ring = true;
    }
    if (atom->getIsAromatic()) {
      profile.aromatic = true;
    }
    if (isDonor(atom)) {
      ++profile.donors;
    }
    if (isAcceptor(atom)) {
      ++profile.acceptors;
    }
  }
  if (profile.volume > 0) {
    profile.hydrophobicity = static_cast<double>(carbonCount) /
                             static_cast<double>(profile.volume);
  }
  profile.amide = containsAmideLike(mol, indices);
  return profile;
}

std::vector<FragmentProfile> buildFragmentProfiles(const ROMol &mol,
                                                   const FeatTreeGraph &tree) {
  std::vector<FragmentProfile> profiles;
  profiles.reserve(boost::num_vertices(tree));
  const auto nodeMap = boost::get(FeatTreeNode_t(), tree);
  for (auto vp = boost::vertices(tree); vp.first != vp.second; ++vp.first) {
    const auto v = *vp.first;
    const UINT_SET &atomSet = nodeMap[v];
    std::set<unsigned int> indices(atomSet.begin(), atomSet.end());
    profiles.push_back(computeFragmentProfile(mol, indices));
  }
  return profiles;
}

std::vector<FragmentProfile> buildFragmentProfiles(const ROMol &mol) {
  auto tree = molToBaseTree(mol);
  if (!tree) {
    return {};
  }
  return buildFragmentProfiles(mol, *tree);
}

double compareFragmentProfiles(const FragmentProfile &lhs,
                               const FragmentProfile &rhs) {
  double score = 0.0;
  int count = 0;

  if (lhs.volume > 0 || rhs.volume > 0) {
    const double diff = std::abs(lhs.volume - rhs.volume);
    const double maxv = static_cast<double>(std::max(lhs.volume, rhs.volume));
    score += 1.0 - diff / maxv;
  } else {
    score += 1.0;
  }
  ++count;

  score += (lhs.ring == rhs.ring) ? 1.0 : 0.0;
  ++count;

  if (lhs.donors > 0 || rhs.donors > 0) {
    const double diff = std::abs(lhs.donors - rhs.donors);
    const double maxv = static_cast<double>(std::max(lhs.donors, rhs.donors));
    score += 1.0 - diff / maxv;
  } else {
    score += 1.0;
  }
  ++count;

  if (lhs.acceptors > 0 || rhs.acceptors > 0) {
    const double diff = std::abs(lhs.acceptors - rhs.acceptors);
    const double maxv =
        static_cast<double>(std::max(lhs.acceptors, rhs.acceptors));
    score += 1.0 - diff / maxv;
  } else {
    score += 1.0;
  }
  ++count;

  score += (lhs.amide == rhs.amide) ? 1.0 : 0.0;
  ++count;

  score += (lhs.aromatic == rhs.aromatic) ? 1.0 : 0.0;
  ++count;

  const double hydDiff = std::abs(lhs.hydrophobicity - rhs.hydrophobicity);
  score += 1.0 - hydDiff;
  ++count;

  return (count > 0) ? score / static_cast<double>(count) : 0.0;
}

namespace {

double directionalSimilarity(const std::vector<FragmentProfile> &lhs,
                             const std::vector<FragmentProfile> &rhs) {
  if (lhs.empty()) {
    return 0.0;
  }
  double sum = 0.0;
  for (const auto &profile : lhs) {
    double best = 0.0;
    for (const auto &candidate : rhs) {
      const double sim = compareFragmentProfiles(profile, candidate);
      if (sim > best) {
        best = sim;
      }
    }
    sum += best;
  }
  return sum / static_cast<double>(lhs.size());
}

}  // namespace

double compareProfileSets(const std::vector<FragmentProfile> &lhs,
                          const std::vector<FragmentProfile> &rhs) {
  if (lhs.empty() || rhs.empty()) {
    return 0.0;
  }
  const double forward = directionalSimilarity(lhs, rhs);
  const double backward = directionalSimilarity(rhs, lhs);
  return 0.5 * (forward + backward);
}

double compareMolecules(const ROMol &lhs, const ROMol &rhs) {
  std::unique_ptr<RWMol> left(new RWMol(lhs));
  std::unique_ptr<RWMol> right(new RWMol(rhs));
  MolOps::sanitizeMol(*left);
  MolOps::sanitizeMol(*right);
  auto leftProfiles = buildFragmentProfiles(*left);
  auto rightProfiles = buildFragmentProfiles(*right);
  return compareProfileSets(leftProfiles, rightProfiles);
}

}  // namespace FeatTrees
}  // namespace RDKit


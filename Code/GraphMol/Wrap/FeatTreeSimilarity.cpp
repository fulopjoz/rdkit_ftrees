//
//  Copyright (C) 2024  RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDBoost/python.h>

#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/RDKitBase.h>

#include <GraphMol/Basement/FeatTrees/FeatTreeSimilarity.h>

#include <memory>

namespace python = boost::python;

namespace {

RDKit::ROMol *copyAndMaybeSanitize(const RDKit::ROMol &mol, bool sanitize) {
  auto *result = new RDKit::ROMol(mol);
  if (sanitize) {
    RDKit::MolOps::sanitizeMol(*result);
  }
  return result;
}

python::dict profileToDict(const RDKit::FeatTrees::FragmentProfile &profile) {
  python::dict data;
  data["volume"] = profile.volume;
  data["ring"] = profile.ring;
  data["donors"] = profile.donors;
  data["acceptors"] = profile.acceptors;
  data["amide"] = profile.amide;
  data["aromatic"] = profile.aromatic;
  data["hydrophobicity"] = profile.hydrophobicity;
  return data;
}

python::list buildProfiles(const RDKit::ROMol &mol, bool sanitize) {
  std::unique_ptr<RDKit::ROMol> work(copyAndMaybeSanitize(mol, sanitize));
  auto profiles = RDKit::FeatTrees::buildFragmentProfiles(*work);
  python::list result;
  for (const auto &profile : profiles) {
    result.append(profileToDict(profile));
  }
  return result;
}

double compareMolecules(const RDKit::ROMol &lhs, const RDKit::ROMol &rhs,
                        bool sanitize) {
  std::unique_ptr<RDKit::ROMol> left(copyAndMaybeSanitize(lhs, sanitize));
  std::unique_ptr<RDKit::ROMol> right(copyAndMaybeSanitize(rhs, sanitize));
  auto leftProfiles = RDKit::FeatTrees::buildFragmentProfiles(*left);
  auto rightProfiles = RDKit::FeatTrees::buildFragmentProfiles(*right);
  return RDKit::FeatTrees::compareProfileSets(leftProfiles, rightProfiles);
}

}  // namespace

BOOST_PYTHON_MODULE(rdftrees) {
  python::scope().attr("__doc__") =
      "Feature-tree utilities for coarse-grained similarity scoring.";

  python::def("build_fragment_profiles", buildProfiles,
              (python::args("mol"), python::args("sanitize") = true),
              "Compute fragment-level feature profiles from the molecule's base feature tree.\n\n"
              "Args:\n"
              "    mol: An RDKit molecule instance.\n"
              "    sanitize: If True (default), a sanitized copy of the molecule is used.\n\n"
              "Returns:\n"
              "    list[dict]: One dictionary per fragment describing volume, donors,\n"
              "    acceptors, aromaticity, amide-like character, and hydrophobicity.");

  python::def("compare_molecules", compareMolecules,
              (python::args("mol1"), python::args("mol2"),
               python::args("sanitize") = true),
              "Compute the symmetric feature-tree similarity between two molecules.\n\n"
              "Args:\n"
              "    mol1: First RDKit molecule.\n"
              "    mol2: Second RDKit molecule.\n"
              "    sanitize: If True (default), sanitized copies of the inputs are used.\n\n"
              "Returns:\n"
              "    float: Similarity score in the range [0, 1].");
}


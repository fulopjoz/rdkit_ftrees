//
//  Copyright (C) 2024  RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/RDLog.h>

#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <GraphMol/Basement/FeatTrees/FeatTreeSimilarity.h>

#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

int main(int argc, char **argv) {
  if (argc < 3) {
    std::cerr << "Usage: ftrees_full <smiles1> <smiles2> [<smiles3> ...]" << std::endl;
    std::cerr << "Computes pairwise feature-tree similarities between the provided SMILES." << std::endl;
    return 1;
  }

  RDLog::InitLogs();

  std::vector<std::unique_ptr<RDKit::RWMol>> molecules;
  molecules.reserve(static_cast<size_t>(argc - 1));
  std::vector<std::string> inputs;
  inputs.reserve(static_cast<size_t>(argc - 1));

  for (int i = 1; i < argc; ++i) {
    const std::string smi = argv[i];
    std::unique_ptr<RDKit::ROMol> parsed(RDKit::SmilesToMol(smi));
    if (!parsed) {
      std::cerr << "Error: failed to parse SMILES '" << smi << "'." << std::endl;
      return 1;
    }
    std::unique_ptr<RDKit::RWMol> mol(new RDKit::RWMol(*parsed));
    try {
      RDKit::MolOps::sanitizeMol(*mol);
    } catch (const RDKit::MolSanitizeException &ex) {
      std::cerr << "Error: sanitization failed for SMILES '" << smi
                << "': " << ex.what() << std::endl;
      return 1;
    }
    molecules.push_back(std::move(mol));
    inputs.push_back(smi);
  }

  std::vector<std::vector<RDKit::FeatTrees::FragmentProfile>> treeProfiles;
  treeProfiles.reserve(molecules.size());

  for (const auto &molPtr : molecules) {
    auto profiles = RDKit::FeatTrees::buildFragmentProfiles(*molPtr);
    if (profiles.empty()) {
      std::cerr << "Error: failed to generate feature tree." << std::endl;
      return 1;
    }
    treeProfiles.push_back(std::move(profiles));
  }

  const size_t count = molecules.size();
  std::cout << std::fixed << std::setprecision(3);
  for (size_t i = 0; i < count; ++i) {
    for (size_t j = i + 1; j < count; ++j) {
      const double similarity =
          RDKit::FeatTrees::compareProfileSets(treeProfiles[i], treeProfiles[j]);
      std::cout << "Similarity(" << inputs[i] << ", " << inputs[j]
                << ") = " << similarity << std::endl;
    }
  }

  return 0;
}


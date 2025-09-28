//
//  Copyright (C) 2024  RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#define NO_IMPORT_ARRAY
#include <RDBoost/python.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/Basement/FeatTrees/FeatTree.h>

#include <sstream>
#include <string>

namespace python = boost::python;

namespace RDKit {
namespace FeatTrees {
namespace {

FeatTreeParams paramsFromObject(const python::object &obj) {
  if (!obj || obj.is_none()) {
    return FeatTreeParams();
  }
  python::extract<FeatTreeParams> extractParams(obj);
  if (extractParams.check()) {
    return extractParams();
  }
  FeatTreeParams params;
  python::extract<python::dict> extractDict(obj);
  if (extractDict.check()) {
    auto dict = extractDict();
    for (python::stl_input_iterator<python::tuple> it(dict.items());
         it != python::stl_input_iterator<python::tuple>(); ++it) {
      const std::string key = python::extract<std::string>((*it)[0]);
      python::object value = (*it)[1];
      if (key == "compressPaths") {
        params.compressPaths = python::extract<bool>(value);
      } else if (key == "mergeBranchGroups") {
        params.mergeBranchGroups = python::extract<bool>(value);
      } else if (key == "maxBranchGroupSize") {
        params.maxBranchGroupSize = python::extract<unsigned int>(value);
      } else if (key == "annotateFeatures") {
        params.annotateFeatures = python::extract<bool>(value);
      } else if (key == "includeZeroNodes") {
        params.includeZeroNodes = python::extract<bool>(value);
      } else if (key == "canonicalize") {
        params.canonicalize = python::extract<bool>(value);
      } else if (key == "ringWeight") {
        params.ringWeight = python::extract<double>(value);
      } else if (key == "connectorWeight") {
        params.connectorWeight = python::extract<double>(value);
      } else if (key == "featureGroupWeight") {
        params.featureGroupWeight = python::extract<double>(value);
      } else {
        throw_value_error(("Unknown FeatTreeParams field: " + key).c_str());
      }
    }
    return params;
  }
  throw_value_error("FeatTreeParams require a mapping or another FeatTreeParams instance");
  return FeatTreeParams();
}

struct PyFeatTree {
  explicit PyFeatTree(FeatTreeGraphSPtr tree) : d_tree(std::move(tree)) {}
  python::list getNodes() const {
    python::list res;
    if (!d_tree) {
      return res;
    }
    const auto nodeMap = boost::get(FeatTreeNode_t(), *d_tree);
    for (auto vp = boost::vertices(*d_tree); vp.first != vp.second; ++vp.first) {
      const auto &node = nodeMap[*vp.first];
      python::dict record;
      python::list atoms;
      for (auto idx : node.atoms) {
        atoms.append(idx);
      }
      record["atoms"] = atoms;
      record["kind"] = static_cast<unsigned int>(node.kind);
      record["flags"] = static_cast<unsigned int>(node.flags);
      record["minRingSize"] = static_cast<unsigned int>(node.minRingSize);
      record["maxRingSize"] = static_cast<unsigned int>(node.maxRingSize);
      record["aromaticAtomCount"] = static_cast<unsigned int>(node.aromaticAtomCount);
      record["heteroAtomCount"] = static_cast<unsigned int>(node.heteroAtomCount);
      res.append(record);
    }
    return res;
  }
  python::list getEdges() const {
    python::list res;
    if (!d_tree) {
      return res;
    }
    const auto edgeMap = boost::get(FeatTreeEdge_t(), *d_tree);
    for (auto ep = boost::edges(*d_tree); ep.first != ep.second; ++ep.first) {
      python::dict data;
      const auto &edge = edgeMap[*ep.first];
      data["ringEndCount"] = static_cast<unsigned int>(edge.ringEndCount);
      data["flags"] = static_cast<unsigned int>(edge.flags);
      res.append(python::make_tuple(static_cast<unsigned int>(
                                    boost::source(*ep.first, *d_tree)),
                                    static_cast<unsigned int>(
                                        boost::target(*ep.first, *d_tree)),
                                    data));
    }
    return res;
  }
  std::string toJSON() const {
    if (!d_tree) {
      return "{}";
    }
    return featTreeToJSON(*d_tree);
  }
  std::string repr() const {
    std::ostringstream oss;
    oss << "<FeatTree nodes=" << boost::num_vertices(*d_tree)
        << " edges=" << boost::num_edges(*d_tree) << ">";
    return oss.str();
  }

  FeatTreeGraphSPtr d_tree;
};

PyFeatTree molToFeatTreePy(const ROMol &mol, const python::object &paramsObj,
                           bool asBaseTree) {
  const auto params = paramsFromObject(paramsObj);
  if (asBaseTree) {
    return PyFeatTree(molToBaseTree(mol, params));
  }
  return PyFeatTree(molToFeatTree(mol, params));
}

double calcFeatTreeSimilarityPy(const ROMol &mol1, const ROMol &mol2,
                                const python::object &paramsObj) {
  const auto params = paramsFromObject(paramsObj);
  return calcFeatTreeSimilarity(mol1, mol2, params);
}

double calcFeatTreeSimilarityGraphPy(const PyFeatTree &tree1,
                                     const PyFeatTree &tree2,
                                     const python::object &paramsObj) {
  const auto params = paramsFromObject(paramsObj);
  PRECONDITION(tree1.d_tree, "Invalid tree");
  PRECONDITION(tree2.d_tree, "Invalid tree");
  return calcFeatTreeSimilarity(*tree1.d_tree, *tree2.d_tree, params);
}

PyFeatTree baseTreeToFeatTreePy(const PyFeatTree &baseTree,
                                const python::object &paramsObj,
                                const ROMol &mol) {
  const auto params = paramsFromObject(paramsObj);
  PRECONDITION(baseTree.d_tree, "Invalid base tree");
  auto copy = FeatTreeGraphSPtr(new FeatTreeGraph(*baseTree.d_tree));
  baseTreeToFeatTree(*copy, params, &mol);
  return PyFeatTree(copy);
}

}  // namespace

void wrapFeatTrees() {
  python::enum_<FeatTreeNodeKind>("FeatTreeNodeKind")
      .value("RingSystem", FeatTreeNodeKind::RingSystem)
      .value("FusedRingSystem", FeatTreeNodeKind::FusedRingSystem)
      .value("Connector", FeatTreeNodeKind::Connector)
      .value("BranchGroup", FeatTreeNodeKind::BranchGroup)
      .value("ZeroNode", FeatTreeNodeKind::ZeroNode)
      .value("FeatureGroup", FeatTreeNodeKind::FeatureGroup);

  python::class_<FeatTreeParams>("FeatTreeParams", "Feature tree construction options",
                                 python::init<>())
      .def_readwrite("compressPaths", &FeatTreeParams::compressPaths)
      .def_readwrite("mergeBranchGroups", &FeatTreeParams::mergeBranchGroups)
      .def_readwrite("maxBranchGroupSize", &FeatTreeParams::maxBranchGroupSize)
      .def_readwrite("annotateFeatures", &FeatTreeParams::annotateFeatures)
      .def_readwrite("includeZeroNodes", &FeatTreeParams::includeZeroNodes)
      .def_readwrite("canonicalize", &FeatTreeParams::canonicalize)
      .def_readwrite("ringWeight", &FeatTreeParams::ringWeight)
      .def_readwrite("connectorWeight", &FeatTreeParams::connectorWeight)
      .def_readwrite("featureGroupWeight", &FeatTreeParams::featureGroupWeight)
      .def("__repr__", [](const FeatTreeParams &params) {
        std::ostringstream oss;
        oss << "FeatTreeParams(compressPaths=" << params.compressPaths
            << ", mergeBranchGroups=" << params.mergeBranchGroups
            << ", maxBranchGroupSize=" << params.maxBranchGroupSize
            << ", annotateFeatures=" << params.annotateFeatures
            << ", includeZeroNodes=" << params.includeZeroNodes
            << ", canonicalize=" << params.canonicalize
            << ", ringWeight=" << params.ringWeight
            << ", connectorWeight=" << params.connectorWeight
            << ", featureGroupWeight=" << params.featureGroupWeight << ')';
        return oss.str();
      });

  python::class_<PyFeatTree>("FeatTree", "Feature tree graph representation",
                             python::init<FeatTreeGraphSPtr>())
      .def("GetNodes", &PyFeatTree::getNodes,
           "Return a list of node dictionaries.")
      .def("GetEdges", &PyFeatTree::getEdges,
           "Return a list of (u, v, data) edge tuples.")
      .def("ToJSON", &PyFeatTree::toJSON, "Serialise the tree to JSON.")
      .def("__repr__", &PyFeatTree::repr);

  python::def("MolToFeatTree", molToFeatTreePy,
              (python::arg("mol"), python::arg("params") = python::object(),
               python::arg("as_base_tree") = false),
              "Constructs a feature tree from a molecule.");
  python::def("BaseTreeToFeatTree", baseTreeToFeatTreePy,
              (python::arg("base_tree"), python::arg("params") = python::object(),
               python::arg("mol")),
              "Applies post-processing stages to a base tree.");
  python::def("CalcFeatTreeSimilarity", calcFeatTreeSimilarityPy,
              (python::arg("mol1"), python::arg("mol2"),
               python::arg("params") = python::object()),
              "Calculates the weighted Jaccard similarity between molecules.");
  python::def("CalcFeatTreeSimilarity", calcFeatTreeSimilarityGraphPy,
              (python::arg("tree1"), python::arg("tree2"),
               python::arg("params") = python::object()),
              "Calculates the weighted Jaccard similarity between feature trees.");
}

}  // namespace FeatTrees
}  // namespace RDKit

BOOST_PYTHON_MODULE(rdFeatTrees) { RDKit::FeatTrees::wrapFeatTrees(); }

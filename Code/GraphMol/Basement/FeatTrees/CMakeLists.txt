if(NOT RDK_BUILD_FEATTREES)
  return()
endif()

set(FeatTrees_sources
    FeatTree.cpp
    FeatTreeUtils.cpp
    FeatTreeSimilarity.cpp)

rdkit_library(FeatTrees
    ${FeatTrees_sources}
    SHARED
    LINK_LIBRARIES
      GraphMol
      MolChemicalFeatures
      RDGeneral)

target_compile_definitions(FeatTrees PRIVATE RDKIT_FEATTREES_BUILD)

rdkit_headers(FeatTree.h
              FeatTreeUtils.h
              FeatTreeSimilarity.h
              DEST GraphMol/Basement/FeatTrees)

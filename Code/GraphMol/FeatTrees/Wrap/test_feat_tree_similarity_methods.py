from rdkit import Chem
from rdkit.Chem import rdFeatTrees


def test_default_signature_compatibility():
    mol = Chem.MolFromSmiles("c1ccccc1")
    params = rdFeatTrees.FeatTreeParams()
    default = rdFeatTrees.CalcFeatTreeSimilarity(mol, mol, params)
    explicit = rdFeatTrees.CalcFeatTreeSimilarity(
        mol, mol, rdFeatTrees.FeatTreeSimilarityMethod.WeightedJaccard, params
    )
    assert abs(default - explicit) < 1e-8


def test_similarity_method_auto_matches_approx():
    mol1 = Chem.MolFromSmiles("CC")
    mol2 = Chem.MolFromSmiles("CCO")
    params = rdFeatTrees.FeatTreeParams()
    params.similarityAutoThreshold = 64
    approx = rdFeatTrees.CalcFeatTreeSimilarity(
        mol1, mol2, rdFeatTrees.FeatTreeSimilarityMethod.ApproxEdit, params
    )
    auto = rdFeatTrees.CalcFeatTreeSimilarity(
        mol1, mol2, rdFeatTrees.FeatTreeSimilarityMethod.Auto, params
    )
    assert abs(approx - auto) < 1e-6


def test_substructure_monotonicity():
    benzene = Chem.MolFromSmiles("c1ccccc1")
    toluene = Chem.MolFromSmiles("Cc1ccccc1")
    butane = Chem.MolFromSmiles("CCCC")
    params = rdFeatTrees.FeatTreeParams()
    sim_bt = rdFeatTrees.CalcFeatTreeSimilarity(
        benzene, toluene, rdFeatTrees.FeatTreeSimilarityMethod.Auto, params
    )
    sim_bb = rdFeatTrees.CalcFeatTreeSimilarity(
        benzene, butane, rdFeatTrees.FeatTreeSimilarityMethod.Auto, params
    )
    sim_td = rdFeatTrees.CalcFeatTreeSimilarity(
        toluene, butane, rdFeatTrees.FeatTreeSimilarityMethod.Auto, params
    )
    assert sim_bt >= sim_bb >= sim_td


def test_identical_molecules_similarity():
    mol = Chem.MolFromSmiles("c1ccccc1")
    params = rdFeatTrees.FeatTreeParams()
    weighted = rdFeatTrees.CalcFeatTreeSimilarity(
        mol, mol, rdFeatTrees.FeatTreeSimilarityMethod.WeightedJaccard, params
    )
    approx = rdFeatTrees.CalcFeatTreeSimilarity(
        mol, mol, rdFeatTrees.FeatTreeSimilarityMethod.ApproxEdit, params
    )
    assert abs(weighted - 1.0) < 1e-8
    assert abs(approx - 1.0) < 1e-8

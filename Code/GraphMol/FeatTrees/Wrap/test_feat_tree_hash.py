from rdkit import Chem
from rdkit.Chem import rdFeatTrees
import json


def test_hash_consistency():
    mol = Chem.MolFromSmiles("c1ccccc1")
    params = rdFeatTrees.FeatTreeParams()
    tree1 = rdFeatTrees.MolToFeatTree(mol, params)
    tree2 = rdFeatTrees.MolToFeatTree(mol, params)
    assert rdFeatTrees.HashFeatTree(tree1) == rdFeatTrees.HashFeatTree(tree2)


def test_hash_consistency_across_pipeline():
    mol = Chem.MolFromSmiles("Cc1ccccc1")
    params = rdFeatTrees.FeatTreeParams()
    base = rdFeatTrees.MolToFeatTree(mol, params, as_base_tree=True)
    full = rdFeatTrees.BaseTreeToFeatTree(base, params, mol)
    assert rdFeatTrees.HashFeatTree(full) == rdFeatTrees.HashFeatTree(base)


def test_json_schema_version_and_params():
    mol = Chem.MolFromSmiles("CCO")
    params = rdFeatTrees.FeatTreeParams()
    params.mergeBranchGroups = False
    tree = rdFeatTrees.MolToFeatTree(mol, params)
    data = json.loads(tree.ToJSON(params))
    assert data["schema_version"] == rdFeatTrees.FEATTREE_SCHEMA_VERSION
    assert data["params"]["similarityAutoThreshold"] == params.similarityAutoThreshold
    assert data["params"]["ringWeight"] == params.ringWeight

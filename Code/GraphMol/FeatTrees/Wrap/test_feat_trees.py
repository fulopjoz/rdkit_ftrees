from rdkit import Chem
from rdkit import Chem
from rdkit.Chem import rdFeatTrees
import json


def _count_kinds(nodes, kind):
    return sum(1 for node in nodes if node["kind"] == int(kind))


def test_mol_to_feat_tree_basic():
    mol = Chem.MolFromSmiles("CCO")
    params = rdFeatTrees.FeatTreeParams()
    tree = rdFeatTrees.MolToFeatTree(mol, params)
    nodes = tree.GetNodes()
    assert nodes
    for node in nodes:
        atoms = node["atoms"]
        assert atoms == sorted(atoms)
    data = json.loads(tree.ToJSON())
    assert len(data["nodes"]) == len(nodes)


def test_feature_annotation_toggle():
    mol = Chem.MolFromSmiles("CCO")
    params = rdFeatTrees.FeatTreeParams()
    params.annotateFeatures = True
    tree = rdFeatTrees.MolToFeatTree(mol, params)
    nodes = tree.GetNodes()
    assert _count_kinds(nodes, rdFeatTrees.FeatTreeNodeKind.FeatureGroup) >= 1


def test_branch_group_toggle():
    mol = Chem.MolFromSmiles("CC(C)C")
    params = rdFeatTrees.FeatTreeParams()
    params.mergeBranchGroups = False
    tree_no_group = rdFeatTrees.MolToFeatTree(mol, params)
    params.mergeBranchGroups = True
    tree_group = rdFeatTrees.MolToFeatTree(mol, params)
    nodes_no_group = tree_no_group.GetNodes()
    nodes_group = tree_group.GetNodes()
    assert _count_kinds(nodes_group, rdFeatTrees.FeatTreeNodeKind.BranchGroup) >= 1
    assert len(nodes_group) <= len(nodes_no_group)


def test_similarity_api():
    mol1 = Chem.MolFromSmiles("c1ccccc1")
    mol2 = Chem.MolFromSmiles("c1ccncc1")
    params = rdFeatTrees.FeatTreeParams()
    sim = rdFeatTrees.CalcFeatTreeSimilarity(mol1, mol2, params)
    assert 0.0 <= sim <= 1.0
    tree1 = rdFeatTrees.MolToFeatTree(mol1, params)
    tree2 = rdFeatTrees.MolToFeatTree(mol2, params)
    sim_graph = rdFeatTrees.CalcFeatTreeSimilarity(tree1, tree2, params)
    assert abs(sim - sim_graph) < 1e-6


def test_base_tree_transformation():
    mol = Chem.MolFromSmiles("CCCCCC")
    params = rdFeatTrees.FeatTreeParams()
    base_tree = rdFeatTrees.MolToFeatTree(mol, params, as_base_tree=True)
    full_tree = rdFeatTrees.BaseTreeToFeatTree(base_tree, params, mol)
    assert len(full_tree.GetNodes()) <= len(base_tree.GetNodes())

import unittest
from upgma import UPGMA_treeConstructor
from Bio import Phylo
import networkx as nx

class TestUPGMA_treeConstructor(unittest.TestCase):

	def test_unsquare(self):
		with self.assertRaises(ValueError):
			UPGMA_treeConstructor('tests/test_unsquared.txt')

	def test_inconsistent_matrix(self):
		with self.assertRaises(ValueError):
			UPGMA_treeConstructor('tests/test_different_length.txt')

	def test_trivial(self):
		tree = UPGMA_treeConstructor('tests/test_trivial.txt').create_tree()

		self.assertEqual(len(tree.clade), 0)

	def test_empty(self):
		tree = UPGMA_treeConstructor('tests/test_empty.txt').create_tree()

		self.assertEqual(tree, None)

	def test_correct(self):
		tree1 = UPGMA_treeConstructor('tests/test.txt')
		tree2 = Phylo.TreeConstruction.DistanceTreeConstructor().upgma(tree1.distances)

		self.assertTrue(nx.is_isomorphic(Phylo.to_networkx(tree1.create_tree()), 
			Phylo.to_networkx(tree2)))


if __name__ == '__main__':
    unittest.main(exit=False)

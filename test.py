import unittest
from upgma import UPGMA_treeConstructor
from Bio.Phylo import TreeConstruction
from Bio.Phylo import Newick

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

	def test_nonexistence(self):
		with self.assertRaises(FileNotFoundError):
			tree = UPGMA_treeConstructor('tests/test_nonexistence.txt').create_tree()

	def test_correct(self):
		tree1 = UPGMA_treeConstructor('test.txt')
		tree2 = TreeConstruction.DistanceTreeConstructor().upgma(tree1.distances)
		print(tree2)
		print(Newick.Tree(tree1.create_tree()))

		self.assertEqual(tree1.create_tree(), tree2)


if __name__ == '__main__':
    unittest.main(exit=False)

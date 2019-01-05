from Bio.Phylo import TreeConstruction
import matplotlib.pyplot as plt
from Bio.Phylo import BaseTree 
from Bio.Phylo import draw as phylo_draw
import copy
#from Bio.Phylo import Newick as BaseTree

def _read_matrix(filename):
	with open(filename) as f:
		header = f.readline().split()

		buf = list(map(lambda l: list(map(float, l[1:].split())), f.readlines()))

	if len(buf) != len(header):
		raise ValueError('Distance matrix must be square')
	
	for i, l in enumerate(buf):
		if len(l) != len(header):
			raise ValueError('Distance matrix must be square')

		buf[i] =  l[i::][::-1]

	if not buf:
		return None

	res = TreeConstruction.DistanceMatrix(names=header[::-1],
			matrix=buf[::-1])
	
	return res

def _join_clades(x, y, dist):
	for i in dist.names:
		if i != x:
			dist[x, i] = (dist[x, i] + dist[y, i]) / 2
		
	dist.names[dist.names.index(x)] = x + y
	del dist[y]

	return dist 

def _find_min(m):
	min_val, i_min, j_min = float('Inf'), 0, 0

	for i in m.names:
		for j in m.names:
			if m[i, j] <= min_val and m[i, j] > 0:
				min_val, i_min, j_min = m[i, j], i, j

	return min_val, i_min, j_min

def _calc_height(clade, height):
	if not clade.clades:
		return height + clade.branch_length

	return height + max([_calc_height(c, clade.branch_length) for c in clade.clades]) 


def _recalc_height(clade, dist):
	if not clade.clades:
		clade.branch_length = dist / 2
	else:
		clade.branch_length = dist / 2 - _calc_height(clade, 0)


class UPGMA_treeConstructor:
	def __init__(self, filename):
		"""
		Creates class instance

		Arguments:
		filename -- the path to file containing
		distance matrix
		"""

		self.distances = _read_matrix(filename)

	def create_tree(self):

		"""Methods that constructs upgma tree
		based on the distance matrix
		"""

		if hasattr(self, 'tree'):
			return self.tree

		if not self.distances:
			self.tree = None
			return None

		clades = [BaseTree.Clade(None, n) for n in self.distances.names]

		find_clade = lambda name: [i for i, el in enumerate(clades) if el.name == name][0]

		while len(self.distances.names) > 1: 
			dist, i, j = _find_min(self.distances)
			i_clade, j_clade = find_clade(i), find_clade(j)
			new_clade = BaseTree.Clade(0, str(i) + str(j))

			_recalc_height(clades[i_clade], dist)
			_recalc_height(clades[j_clade], dist)

			new_clade.clades.append(clades[i_clade])
			new_clade.clades.append(clades[j_clade])

			if j_clade > i_clade:
				i_clade, j_clade = j_clade, i_clade

			clades.pop(i_clade)
			clades.pop(j_clade)

			clades.append(new_clade)
			self.distances = _join_clades(i, j, self.distances)

		self.tree = BaseTree.Tree(clades[0])

		return self.tree

	def draw(self,
		filename=None):

		"""Method that draws the tree to the
		file with given name or shows it in 
		a pop-up window

		Keyword arguments:
		filename -- the name of file to save the result to
		"""

		if not hasattr(self, 'tree'):
			self.create_tree()

		show = True if not filename else False
		phylo_draw(self.tree,
			do_show=show)

		if filename:
			plt.savefig(filename)
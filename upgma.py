from Bio.Phylo import TreeConstruction
import matplotlib.pyplot as plt

def _read_matrix(filename):
	with open(filename) as f:
		header = f.readline().split()

		buf = list(map(lambda l: list(map(float, l[1:].split())), f.readlines()))

	if len(buf) != len(header):
		raise ValueError('Distance matrix must be square')
	
	#TODO prettify this block
	for i, l in enumerate(buf):
		if len(l) != len(header):
			raise ValueError('Distance matrix must be square')

		buf[i] =  l[i::]
		buf[i] = buf[i][::-1]


	if not buf:
		return None

	res = TreeConstruction.DistanceMatrix(names=header[::-1],
			matrix=buf[::-1])
	
	return res

def _join_clades(x, y, sizex, sizey, dist):
	
	for i in dist.names:
		dist[i, x] = (sizex * dist[i, x] + sizey * dist[i, y]) / (sizex + sizey)
		
	dist.names[dist.names.index(x)] = x + y
	del dist[y]

	return dist 

def _find_min(m):
	min_val, i_min, j_min = float('Inf'), 0, 0

	for i in m.names:
		for j in m.names:
			if m[i, j] < min_val and m[i, j] > 0:
				min_val, i_min, j_min = m[i, j], i, j

	return min_val, i_min, j_min 


def _calc_height(clade, dist):
	if len(clade.clades) == 0:
		clade.branch_length = dist / 2
	else:
		clade.branch_length = dist / 2 - clade.clades[-1].branch_length


class UPGMA_treeConstructor:
	def __init__(self, filename):
		self.distances = _read_matrix(filename)

	def create_tree(self):
		if not self.distances:
			self.tree = None
			return None

		clades = [Phylo.BaseTree.Clade(None, n) for n in self.distances.names]

		find_clade = lambda name: [i for i, el in enumerate(clades) if el.name == name][0]

		while len(self.distances.names) > 1: 
			dist, i, j = _find_min(self.distances)
			i_clade, j_clade = find_clade(i), find_clade(j)
			new_clade = Phylo.BaseTree.Clade(None, str(i) + str(j))

			_calc_height(clades[i_clade], dist)
			_calc_height(clades[j_clade], dist)

			new_clade.clades.append(clades[i_clade])
			new_clade.clades.append(clades[j_clade])

			if j_clade > i_clade:
				i_clade, j_clade = j_clade, i_clade

			clades.pop(i_clade)
			clades.pop(j_clade)

			clades.append(new_clade)
			self.distances = _join_clades(i, j, len(new_clade.clades), len(new_clade.clades), self.distances)

		self.tree = Phylo.BaseTree.Tree(clades[0])

		return self.tree

	def draw(self,
		filename=None):

		if not hasattr(self, 'tree'):
			self.create_tree()

		show = True if not filename else False
		Phylo.draw(self.tree,
			do_show=show)

		if filename:
			plt.savefig(filename)




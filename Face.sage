class Face:
	def __init__(self):
		self.Dimension = -1
		self.InnerNormals = []
		self.Vertices = set()
		self.Neighbors = set()
		self.Children = set()
		self.Parents = set()
		self.MyCone = 0
		self.CPolyhedron = 0
		self.Equations = 0
		self.Inequalities = 0
	def PrintProps(self):
		print ""
		print "Below is a face of dimension " , self.Dimension
		print "InnerNormal(s) = " , self.InnerNormals
		print "Vertices = ", self.Vertices
		print "Neighbors = ", self.Neighbors
		print "Children = ", self.Children
		print "Parents = ", self.Parents
		print "Cone = ", self.MyCone
		print "CPolyhedron = ", self.CPolyhedron
		print ""
		return


class ReducedConeClass:
	def __init__(self):
		self.MyCone = Cone([[0]])
		self.ReducedCone = Cone([[0]])
		self.EdgeTuples = []
		self.Reductions = []
		self.DimensionsRemoved = []
		self.CPolyhedron = 0
		self.Has_CPolyhedron = False
	def __repr__(self):
		L = []
		for ray in self.MyCone.rays():
			L.append(ray)
		L.sort()
		LL = self.DimensionsRemoved
		LL.sort()
		return str(L) + str(LL)
	def __hash__(self):
		return hash(self.__repr__())
	def __eq__(self, other):
		if isinstance(other, ReducedConeClass):
			return (self.__repr__() == other.__repr__())
		else:
			return False
	def __ne__(self, other):
		return (not self.__eq__(other))

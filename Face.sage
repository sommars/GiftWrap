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

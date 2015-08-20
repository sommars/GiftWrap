class Face:
	def __init__(self):
		self.InnerNormal = -1
		self.KnowsAllNeighbors = False
		self.Parents = []
		self.Children = []
		self.Neighbors = []
		self.Vertices = []
		self.Dimension = -1
	def PrintProps(self):
		print ""
		print "Below is a face of dimension " , self.Dimension
		print "InnerNormal = " , self.InnerNormal
		print "KnowsAllNeighbors = " , self.KnowsAllNeighbors
		print "Vertices = ", self.Vertices
		print ""
		return

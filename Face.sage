class Face:
	def __init__(self):
		self.InnerNormal = -1
		self.Parents = []
		self.Children = []
		self.Neighbors = []
		self.Vertices = []
		self.Dimension = -1
	def PrintProps(self):
		print ""
		print "Below is a face of dimension " , self.Dimension
		print "InnerNormal = " , self.InnerNormal
		print "Vertices = ", self.Vertices
		print "Dimension = ", self.Dimension
		print "Neighbors = ", self.Neighbors
		print "Parents = ", self.Parents
		print "Children = ", self.Children
		print ""
		return

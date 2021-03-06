load("GiftWrap.sage")

def NormVector(V):
	"""
	Finds the norm of a vector
	"""
	Denom = 0
	for i in xrange(len(V)):
		Denom += V[i]^2
	Denom = Denom^(1/2)
	return [V[i]/Denom for i in xrange(len(V))]

#-------------------------------------------------------------------------------
def InnerProduct(V1,V2):
	"""
	Finds the inner product of two vectors
	"""
	InnerProd = 0
	for i in xrange(len(V1)):
		InnerProd += V1[i]*V2[i]
	return InnerProd

#-------------------------------------------------------------------------------
def FindInitialForm(Pts, V):
	"""
	Computes the initial form of a vector and a set of points.
	"""
	InForm = []
	MinimalIP = 'Test'
	for Pt in Pts:
		IP = InnerProduct(V,Pt)
		if MinimalIP == 'Test':
			MinimalIP = IP
			InForm = [Pt]
		elif MinimalIP > IP:
			MinimalIP = IP
			InForm = [Pt]
		elif IP == MinimalIP:
			InForm.append(Pt)
	InForm.sort()
	return InForm

#-------------------------------------------------------------------------------
def ComputeFacetPretropisms(AFaces, BFaces):
	"""
	Testing function
	"""
	#Currently is an inefficient implementation
	AFacetVectors = [(NormVector(AFaces[1][i].InnerNormals[0]),i) for i in xrange(len(AFaces[1]))]
	BFacetVectors = [NormVector(BFaces[1][i].InnerNormals[0]) for i in xrange(len(BFaces[1]))]

	FacetPretropisms = []
	NotFacetPretropisms = []
	for Vector in AFacetVectors:
		if Vector[0] in BFacetVectors:
			FacetPretropisms.append(Vector)
		else:
			NotFacetPretropisms.append(Vector)
	return FacetPretropisms, NotFacetPretropisms

#-------------------------------------------------------------------------------
def FasterWrap(PolyAsPts):
	"""
	Wraps Sage's convex hull output to find only the edge skeleton where each 
	edge knows its neighbors and knows its cone.
	"""

	def GetCPolyhedron(OldCone):
		cs = Constraint_System()
		BoolList = [True for i in xrange(OldCone.lattice_dim())]
		OldPolyhedron = OldCone.polyhedron()
		for Ineq in [i[1:] for i in OldPolyhedron.inequalities_list()]:
			Expr = 0
			for i in xrange(len(Ineq)):
				Expr = Expr + Ineq[i]*Variable(i)
			cs.insert(Expr >= 0)
		for Eq in [i[1:] for i in OldPolyhedron.equations_list()]:
			Expr = 0
			for i in xrange(len(Eq)):
				Expr = Expr + Eq[i]*Variable(i)
			cs.insert(Expr == 0)
		return C_Polyhedron(cs)

	Faces = [[],[]]
	SagePolyhedron = Polyhedron(PolyAsPts)
	PositiveLineality = SagePolyhedron.equations_list()
	PositiveLineality = [PositiveLineality[i][1:] for i in xrange(len(PositiveLineality))]
	Lineality = PositiveLineality + [[-i for i in j] for j in PositiveLineality]
	
	IndexToPointMap = {}
	PointToIndexMap = {}
	NecessaryVertices = []
	for i in xrange(len(SagePolyhedron.vertices())):
		Vertex = list(SagePolyhedron.vertices()[i])
		IndexToPointMap[i] = tuple(Vertex)
		PointToIndexMap[tuple(Vertex)] = i
		NecessaryVertices.append(Vertex)

	for i in xrange(len(SagePolyhedron.faces(1))):
		NewFace = Face()
		NewFace.Dimension = 1
		Vertices = SagePolyhedron.faces(1)[i].vertices()
		for Vertex in Vertices:
			NewFace.Vertices.add(PointToIndexMap[tuple(Vertex)])
		Faces[0].append(NewFace)

	InitialFormTuples = []
	for Ineq in SagePolyhedron.inequalities_list():
		InitialFormTuples.append((Ineq[1:],FindInitialForm(NecessaryVertices, Ineq[1:])))

	Dim = SagePolyhedron.dim()-1
	for i in xrange(len(SagePolyhedron.faces(Dim))):
		NewFace = Face()
		Vertices = SagePolyhedron.faces(Dim)[i].vertices()
		NewFace.Vertices = set([PointToIndexMap[tuple(Vertex)] for Vertex in Vertices])

		FullDimVertices = [list(Vertex) for Vertex in Vertices]
		FullDimVertices.sort()
		for Tuple in InitialFormTuples:
			if Tuple[1] == FullDimVertices:
				NewFace.InnerNormals.append(Tuple[0])
				break
		Faces[1].append(NewFace)

	for i in xrange(len(Faces[0])):
		Neighbor1 = Faces[0][i]
		Neighbor1.InnerNormals = copy(Lineality)
		for j in xrange(len(Faces[0])):
			if i != j:
				Neighbor2 = Faces[0][j]
				if len(Neighbor1.Vertices.intersection(Neighbor2.Vertices)) == 1:
					Neighbor1.Neighbors.add((j,0))
					Neighbor2.Neighbors.add((i,0))
		for j in xrange(len(Faces[1])):
			PossibleParent = Faces[1][j]
			if len(Neighbor1.Vertices.intersection(PossibleParent.Vertices)) == 2:
				Neighbor1.InnerNormals += PossibleParent.InnerNormals

	for i in xrange(len(Faces[0])):
		Faces[0][i].MyCone = Cone(Faces[0][i].InnerNormals)
		Faces[0][i].CPolyhedron = GetCPolyhedron(Faces[0][i].MyCone)

	print "Lineality space =", Lineality
	return Faces, IndexToPointMap, PointToIndexMap, NecessaryVertices

load("GiftWrap.sage")

def NormVector(V):
	Denom = 0
	for i in xrange(len(V)):
		Denom += V[i]^2
	Denom = Denom^(1/2)
	return [V[i]/Denom for i in xrange(len(V))]

#-------------------------------------------------------------------------------
def InnerProduct(V1,V2):
	InnerProd = 0
	for i in xrange(len(V1)):
		InnerProd += V1[i]*V2[i]
	return InnerProd

#-------------------------------------------------------------------------------
def FindInitialForm(Pts, V):
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
	return InForm

#-------------------------------------------------------------------------------
def ComputeFacetPretropisms(AFaces, BFaces):
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
	Faces = [[],[]]
	SagePolyhedron = Polyhedron(PolyAsPts)
	IndexToPointMap = {}
	PointToIndexMap = {}
	NecessaryVertices = []
	for i in xrange(len(SagePolyhedron.vertices())):
		Vertex = list(SagePolyhedron.vertices()[i])
		IndexToPointMap[i] = tuple(Vertex)
		PointToIndexMap[tuple(Vertex)] = i
		NecessaryVertices.append(Vertex)		

	Barycenter = FindBarycenter(NecessaryVertices)

	for i in xrange(len(SagePolyhedron.faces(1))):
		NewFace = Face()
		NewFace.Dimension = 1
		Vertices = SagePolyhedron.faces(1)[i].vertices()
		for Vertex in Vertices:
			NewFace.Vertices.add(PointToIndexMap[tuple(Vertex)])
		Faces[0].append(NewFace)

	Dim = SagePolyhedron.dim()-1
	for i in xrange(len(SagePolyhedron.faces(Dim))):
		NewFace = Face()
		NewFace.Dimension = Dim
		Vertices = SagePolyhedron.faces(Dim)[i].vertices()
		for Vertex in Vertices:
			NewFace.Vertices.add(PointToIndexMap[tuple(Vertex)])
		FullDimVertices = Vertices
		Normal = GetNormalFromHNF(GetHNF(FullDimVertices))
		NewFace.InnerNormals = [MakeNormalPointInDirectionOfPt(Normal, Barycenter, FullDimVertices[0])]
		Faces[1].append(NewFace)

	for i in xrange(len(Faces[0])):
		Neighbor1 = Faces[0][i]
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

	Startt = time()
	for i in xrange(len(Faces[0])):
		Faces[0][i].MyCone = Cone(Faces[0][i].InnerNormals)
	return Faces, IndexToPointMap, PointToIndexMap, NecessaryVertices

#-------------------------------------------------------------------------------
def WrapHull(PolyAsPts):
	Faces = []
	SagePolyhedron = Polyhedron(PolyAsPts)
	IndexToPointMap = {}
	PointToIndexMap = {}
	NecessaryVertices = []
	for i in xrange(len(SagePolyhedron.vertices())):
		Vertex = list(SagePolyhedron.vertices()[i])
		IndexToPointMap[i] = tuple(Vertex)
		PointToIndexMap[tuple(Vertex)] = i
		NecessaryVertices.append(Vertex)

	Barycenter = FindBarycenter(NecessaryVertices)

	for i in xrange(1,SagePolyhedron.dim()):
		Faces.append([])
		for j in xrange(len(SagePolyhedron.faces(i))):
			NewFace = Face()
			NewFace.Dimension = i
			Vertices = SagePolyhedron.faces(i)[j].vertices()
			for Vertex in Vertices:
				NewFace.Vertices.add(PointToIndexMap[tuple(Vertex)])
			if i == SagePolyhedron.dim() - 1:
				FullDimVertices = []
				for Vertex in NewFace.Vertices:
					FullDimVertices.append(IndexToPointMap[Vertex])
				Normal = GetNormalFromHNF(GetHNF(FullDimVertices))
				NewFace.InnerNormals = [MakeNormalPointInDirectionOfPt(Normal, Barycenter, FullDimVertices[0])]
			Faces[i-1].append(NewFace)

	for i in xrange(0,len(Faces)-1):
		for j in xrange(len(Faces[i])):
			ChildFace = Faces[i][j]
			for k in xrange(len(Faces[i+1])):
				ParentFace = Faces[i+1][k]
				if ChildFace.Vertices.issubset(ParentFace.Vertices):
					ChildFace.Parents.add(k)
					ParentFace.Children.add(j)
	
	for i in xrange(1,len(Faces)):
		Index = len(Faces) - i
		for j in xrange(len(Faces[Index])):
			NewFace = Faces[Index][j]
			for Child in NewFace.Children:
				Child = Faces[Index-1][Child] 
				Child.InnerNormals += NewFace.InnerNormals


	for i in xrange(len(Faces[0])):
		for j in xrange(len(Faces[0])):
			if i != j:
				Neighbor1 = Faces[0][i]
				Neighbor2 = Faces[0][j]
				if len(Neighbor1.Vertices.intersection(Neighbor2.Vertices)) == 1:
					Neighbor1.Neighbors.add((j,0))
					Neighbor2.Neighbors.add((i,0))

	for i in xrange(len(Faces[0])):
		Faces[0][i].MyCone = Cone(Faces[0][i].InnerNormals)

	return Faces, IndexToPointMap, PointToIndexMap, NecessaryVertices

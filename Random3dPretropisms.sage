def DoRandomPretropismTest():
	R.<x,y,z> = QQ[]
	A = R.random_element(10,15)
	B = R.random_element(10,15)
	Rays = R.ideal([A, B]).groebner_fan().tropical_intersection().rays()
	print "Polynomial1 is ", A
	print "Polynomial2 is ", B
	A = [[Integer(j) for j in i] for i in A.exponents()]
	B = [[Integer(j) for j in i] for i in B.exponents()]
	print A
	print B


	AFaces, AIndexToPointMap, APointToIndexMap, A = GiftWrap(A)
	BFaces, BIndexToPointMap, BPointToIndexMap, B = GiftWrap(B)


	FacetPretropisms, NotFacetPretropisms = ComputeFacetPretropisms(AFaces, BFaces)

	if len(FacetPretropisms) > 0:
		print "Polytopes have a common factor, so they are not sufficiently random"
		return

	ConeSet = set()
	EdgePretropisms = []
	for AEdge in AFaces[0]:
		Normal = [AEdge.InnerNormals[0][i] + AEdge.InnerNormals[1][i] for i in xrange(len(AEdge.InnerNormals[0]))]
		BInitialForm = FindInitialForm(B, Normal)
	
		EdgesToTest = set()
		BInitialIndices = set([BPointToIndexMap[tuple(Pt)] for Pt in BInitialForm])
		if len(BInitialIndices) == 1:
			for i in xrange(len(BFaces[0])):
				BFace = BFaces[0][i]
				if BInitialIndices.issubset(BFace.Vertices):
					EdgesToTest.add(i)			
		elif len(BInitialIndices) == 2:
			for i in xrange(len(BFaces[0])):
				BFace = BFaces[0][i]
				if BInitialIndices == BFace.Vertices:
					EdgesToTest.add(i)
					break
		else:
			for i in xrange(len(BFaces[0])):
				BFace = BFaces[0][i]
				if len(BInitialIndices.intersection(BFace.Vertices)) == 2:
					EdgesToTest.add(i)
		PretropEdges = set()
		NotPretropEdges = set()
		while(len(EdgesToTest) > 0):
			TestEdge = EdgesToTest.pop()
			BEdge = BFaces[0][TestEdge]
			if ConesDoIntersect(BEdge.InnerNormals, AEdge.InnerNormals):
				PretropEdges.add(TestEdge)
				for Neighbor in BEdge.Neighbors:
					if Neighbor[0] not in PretropEdges and Neighbor[0] not in NotPretropEdges:
						EdgesToTest.add(Neighbor[0])
			else:
				NotPretropEdges.add(TestEdge)

		#We need to test if all of these edges make up a facet.
		PretropFaces = set()
		for i in xrange(len(BFaces[1])):
			if BFaces[1][i].Children.issubset(EdgesToTest):
				PretropFaces.add(i)
		for Face in PretropFaces:
				EdgesToTest = EdgesToTest.difference(BFaces[1][Face].Children)
		ACone = Cone(AEdge.InnerNormals)
		for Edge in PretropEdges:
			ConeSet.add(ACone.intersection(Cone(BFaces[0][Edge].InnerNormals)).rays())
		for Face in PretropFaces:
			ConeSet.add(ACone.intersection(Cone(BFaces[1][Face].InnerNormals)).rays())


	ConeList = list(ConeSet)
	ConeList.sort()
	for NewCone in ConeList:
		print NewCone
	for NewRay in Rays:
		print NewRay
	return len(ConeList), len(Rays)

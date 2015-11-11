load("Pretropism_Util.sage")
from time import time
def DoRandomPretropismTest(nvars):
	def IntersectCones(Index, NewCone):
		global ConeSet
		global IntersectingRefList
		global HullInfoMap
		Faces = HullInfoMap[(Index + 1,"Faces")]
		for i in IntersectingRefList[Index]:
			TempCone = NewCone.intersection(Cone(Faces[i[0]][i[1]].InnerNormals))
			if TempCone.dim() > 0:
				if Index == len(IntersectingRefList) - 1 and len(TempCone.rays()) == 1:
					ConeSet.add(TempCone.rays())
				else:
					IntersectCones(Index+1,TempCone)
		return
	PolyString = ""
	for i in xrange(nvars - 1):
		PolyString += "x_" + str(i) + ','
	PolyString += "x_" + str(nvars - 1)
	R = PolynomialRing(QQ, nvars, PolyString)
	Polys = [R.random_element(100,10) for i in xrange(nvars-1)]

	for i in xrange(len(Polys)):
		print "Polynomial",  i, " is ", Polys[i]
	starttime = time()
	Rays = R.ideal(Polys).groebner_fan().tropical_intersection().rays()
	GFanTime = time() - starttime

	PolysAsPts = [[[Integer(j) for j in i] for i in Poly.exponents()] for Poly in Polys]

	global HullInfoMap
	HullInfoMap = {}
	HullTime = time()
	for i in xrange(len(PolysAsPts)):
		print PolysAsPts[i]
		Faces, IndexToPointMap, PointToIndexMap, Pts = GiftWrap(PolysAsPts[i])
		HullInfoMap[(i,"Faces")] = Faces
		HullInfoMap[(i,"IndexToPointMap")] = IndexToPointMap
		HullInfoMap[(i,"PointToIndexMap")] = PointToIndexMap
		HullInfoMap[(i,"Pts")] = Pts
	HullTime = time() - HullTime
	JeffStartTime = time()

	global ConeSet
	ConeSet = set()
	EdgePretropisms = []
	AFaces = HullInfoMap[(0,"Faces")]
	AIndexToPointMap =	HullInfoMap[(0,"IndexToPointMap")]
	APointToIndexMap = HullInfoMap[(0,"PointToIndexMap")]
	A =	HullInfoMap[(0,"Pts")]
	FacetPretropisms, NotFacetPretropisms = ComputeFacetPretropisms(AFaces, HullInfoMap[(1,"Faces")])
	if len(FacetPretropisms) > 0:
		print "Polytopes have a common factor, so they are not sufficiently random"
		return

	for AEdge in AFaces[0]:
		global IntersectingRefList
		IntersectingRefList = []
		Normal = [AEdge.InnerNormals[0][i] + AEdge.InnerNormals[1][i] for i in xrange(len(AEdge.InnerNormals[0]))]
		for PolytopeIndex in xrange(1,len(PolysAsPts)):
			BFaces = HullInfoMap[(PolytopeIndex,"Faces")]
			BIndexToPointMap =	HullInfoMap[(PolytopeIndex,"IndexToPointMap")]
			BPointToIndexMap = HullInfoMap[(PolytopeIndex,"PointToIndexMap")]
			B =	HullInfoMap[(PolytopeIndex,"Pts")]
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
			IndexOfIntRefList = len(IntersectingRefList)
			IntersectingRefList.append([])

			for i in xrange(1, len(BFaces)):
				BFacesIndex = len(BFaces) - i
				FacesToKnockOut = []
				for j in xrange(len(BFaces[BFacesIndex])):
					if BFaces[BFacesIndex][j].Children.issubset(EdgesToTest):
						IntersectingRefList[IndexOfIntRefList].add((BFacesIndex, j))
						FacesToKnockOut.append(j)
				for Face in FacesToKnockOut:
					PretropEdges = PretropEdges.difference(BFaces[BFacesIndex][Face].Children)
			for Edge in PretropEdges:
				IntersectingRefList[IndexOfIntRefList].append((0,Edge))

		#What we want here is a list of lists. Each element in the list should look like: (indextodimensiontolookinto, indextowhichelementinBFaces[i][?]

		ACone = Cone(AEdge.InnerNormals)
		IntersectCones(0, ACone)

	ConeList = list(ConeSet)
	ConeList.sort()
	for NewCone in ConeList:
		print NewCone#, len(FindInitialForm(HullInfoMap[(0,"Pts")],NewRay)),len(FindInitialForm(HullInfoMap[(1,"Pts")],NewRay)),len(FindInitialForm(HullInfoMap[(2,"Pts")],NewRay)),len(FindInitialForm(HullInfoMap[(3,"Pts")],NewRay)),len(FindInitialForm(HullInfoMap[(4,"Pts")],NewRay))
	for NewRay in Rays:
		NewRay = [-NewRay[i] for i in xrange(len(NewRay))]
		print NewRay#, len(FindInitialForm(HullInfoMap[(0,"Pts")],NewRay)),len(FindInitialForm(HullInfoMap[(1,"Pts")],NewRay)),len(FindInitialForm(HullInfoMap[(2,"Pts")],NewRay)),len(FindInitialForm(HullInfoMap[(3,"Pts")],NewRay))
	print "GFANTIME", GFanTime
	print "HullTime", HullTime
	print "JeffTime", time() - JeffStartTime
	return len(ConeList), len(Rays)
	
	

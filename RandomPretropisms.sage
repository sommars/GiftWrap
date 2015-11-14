load("Pretropism_Util.sage")
from time import time
#-------------------------------------------------------------------------------
def DoTests(nvars):
	PolyString = ""
	for i in xrange(nvars - 1):
		PolyString += "x_" + str(i) + ','
	PolyString += "x_" + str(nvars - 1)
	R = PolynomialRing(QQ, nvars, PolyString)
	HighestExp = 100
	NumberOfTerms = 10
	Polys = [R.random_element(HighestExp,NumberOfTerms) for i in xrange(nvars-1)]
	for i in xrange(len(Polys)):
		print "Polynomial",  i, " is ", Polys[i]
	PolysAsPts = [[[Integer(j) for j in i] for i in Poly.exponents()] for Poly in Polys]


	global HullInfoMap
	HullInfoMap = {}
	HullTime = time()
	PtsList = []
	for i in xrange(len(PolysAsPts)):
		print PolysAsPts[i]
		Faces, IndexToPointMap, PointToIndexMap, Pts = GiftWrap(PolysAsPts[i],True)
		HullInfoMap[(i,"Faces")] = Faces
		HullInfoMap[(i,"IndexToPointMap")] = IndexToPointMap
		HullInfoMap[(i,"PointToIndexMap")] = PointToIndexMap
		HullInfoMap[(i,"Pts")] = Pts
		PtsList.append(Pts)
		if Faces == 0 or len(Faces) != nvars - 1:
			print "Not all polytopes are the same dimension"
			return 0, 0
	HullTime = time() - HullTime
	for i in xrange(5):
		print ""
	print "ConvexHullTime", HullTime
	print ""

	DoGfan(Polys,R)
	#DoRandomPretropismTest(PolysAsPts, HullInfoMap)
	#DoMinkowskiSum(Polys, PtsList)
	DoCayleyPolytope(PtsList)
	#DoNaiveAlgorithm(HullInfoMap)
	return



#-------------------------------------------------------------------------------
def DoRandomPretropismTest(PolysAsPts, HullInfoMap):
	def IntersectCones(Index, NewCone):
		global ConeSet
		global IntersectingRefList
		global HullInfoMap
		Faces = HullInfoMap[(Index + 1,"Faces")]
		for i in IntersectingRefList[Index]:
			TempCone = NewCone.intersection(Cone(Faces[i[0]][i[1]].InnerNormals))
			if TempCone.dim() > 0:
				if Index == len(IntersectingRefList) - 1:
					if len(TempCone.rays()) == 1:
						for Ray in TempCone.rays():
							ConeSet.add(tuple(Ray.list()))
				else:
					IntersectCones(Index+1,TempCone)
		return
	NewAlgStart = time()

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
		Normal = [0 for i in xrange(len(AEdge.InnerNormals[0]))]
		for InnerNormal in AEdge.InnerNormals:
			Normal = [Normal[i] + InnerNormal[i] for i in xrange(len(Normal))]
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
	global Rays
	if len(ConeList) == len(Rays):
		for i in xrange(len(ConeList)):
			if list(ConeList[i]) != Rays[i]:
				print "UNEQUAL!", Rays[i], ConeList[i]
				return 0, 0
	else:
		print "Lists are unequal lengths", len(ConeList), len(Rays)
		return 0, 0

	print "NewAlg took", time() - NewAlgStart, "seconds."
	print "NewAlg found", len(ConeList), "rays."
	return# len(ConeList), time() - NewAlgStart#, len(HullInfoMap[(0,"Faces")][0]), len(HullInfoMap[(1,"Faces")][0]), len(HullInfoMap[(2,"Faces")][0])
	
#-------------------------------------------------------------------------------
def DoMinkowskiSum(Polys, PtsList):
	MinkowskiStart = time()
	ProdPoly = 1
	for Poly in Polys:
		ProdPoly = ProdPoly*Poly
	ProdPolysAsPts = [[Integer(j) for j in i] for i in ProdPoly.exponents()]

	Faces, IndexToPointMap, PointToIndexMap, Pts = GiftWrap(ProdPolysAsPts,True)
	Normals = []
	for Face in Faces[len(Faces)- 1]:
		Normals.append(Face.InnerNormals[0])
	Counter = 0
	for Normal in Normals:
		ShouldIncCounter = True
		for Pts in PtsList:
			if len(FindInitialForm(Pts, Normal)) < 2:
				ShouldIncCounter = False
				break
		if ShouldIncCounter == True:
			Counter += 1
	print "Minkowski took", time() - MinkowskiStart, "seconds."
	print "Minkowski found", Counter, "rays."
	return

#-------------------------------------------------------------------------------
def DoCayleyPolytope(PtsList):
	CayleyStart = time()
	NewPtSet = []
	
	#This works in 3d, but not in nd
	for i in xrange(len(PtsList)):
		for j in xrange(len(PtsList[i])):
			NewPtSet.append(PtsList[i][j] +[i])
	print NewPtSet
	for i in xrange(len(NewPtSet)):
		for j in xrange(len(NewPtSet[j])):
			NewPtSet[i][j] = Integer(NewPtSet[i][j])

	Faces, IndexToPointMap, PointToIndexMap, Pts = GiftWrap(NewPtSet,True)

	print "PTSLIST BELOW"
	print PtsList[0]
	print PtsList[1]
	NormalList = []
	for Face in Faces[len(Faces) - 1]:
		VertexList = []
		for Pt in Face.Vertices:
			VertexList.append(IndexToPointMap[Pt][:-1])
		Count1 = 0
		Count2 = 0
		for Vertex in VertexList:
			if Vertex in PtsList[0]:
				Count1 += 1
			if Vertex in PtsList[1]:
				Count2 += 1
		if Count1 > 1 and Count2 > 1:
			NormalList.append(Face.InnerNormals[0])
	NormalList.sort()

	print "GFANRAYS"
	for Ray in Rays:
		print Ray
	print ""
	for Normal in NormalList:
		print Normal
	print "Cayley found", len(NormalList), "rays."
	print "CayleyTime", time() - CayleyStart
	return 

#-------------------------------------------------------------------------------
def DoGfan(Polys, R):
	starttime = time()
	global Rays
	Rays = R.ideal(Polys).groebner_fan().tropical_intersection().rays()
	for i in xrange(len(Rays)):
		Rays[i] = [-Rays[i][j] for j in xrange(len(Rays[i]))]
	Rays.sort()
	GfanTime = time() - starttime
	print "Gfan took", GfanTime, "seconds."
	print "Gfan found", len(Rays), "rays."
	return
	
#-------------------------------------------------------------------------------
def DoNaiveAlgorithm(HullInfoMap):
	print "Naive algorithm not yet implemented."
	return

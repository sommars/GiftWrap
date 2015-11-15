load("Pretropism_Util.sage")
from time import time
#-------------------------------------------------------------------------------
def DoTests(nvars):
	PolyString = ""
	for i in xrange(nvars - 1):
		PolyString += "x_" + str(i) + ','
	PolyString += "x_" + str(nvars - 1)
	R = PolynomialRing(QQ, nvars, PolyString)
	HighestExp = 1000
	NumberOfTerms = 15
	Polys = [R.random_element(HighestExp,NumberOfTerms) for i in xrange(nvars-1)]

	for i in xrange(len(Polys)):
		print "Polynomial",  i, " is ", Polys[i]
	PolysAsPts = [[[Integer(j) for j in i] for i in Poly.exponents()] for Poly in Polys]


	global HullInfoMap
	HullInfoMap = {}
	HullTime = time()
	PtsList = []
	EdgeProd = 1
	for i in xrange(len(PolysAsPts)):
		#print PolysAsPts[i]
		Faces, IndexToPointMap, PointToIndexMap, Pts = GiftWrap(PolysAsPts[i],True)
		HullInfoMap[(i,"Faces")] = Faces
		HullInfoMap[(i,"IndexToPointMap")] = IndexToPointMap
		HullInfoMap[(i,"PointToIndexMap")] = PointToIndexMap
		HullInfoMap[(i,"Pts")] = Pts
		EdgeProd = EdgeProd*len(Faces[0])
		PtsList.append(Pts)
		if Faces == 0 or len(Faces) != nvars - 1:
			print "Not all polytopes are the same dimension"
			return 0, 0

	print ""
	HullTime = time() - HullTime
	print "ConvexHullTime", HullTime
	#DoGfan(Polys,R)
	TempStr = DoNewAlgorithm(PolysAsPts, HullInfoMap)
	#DoMinkowskiSum(Polys, PtsList)
	#DoCayleyPolytope(PtsList)
	#DoNaiveAlgorithm(HullInfoMap, nvars - 1)
	return str(nvars) + ',' + str(EdgeProd) +',' + TempStr

#-------------------------------------------------------------------------------
def DoNewAlgorithm(PolysAsPts, HullInfoMap):
	def IntersectCones(Index, NewCone):
		global ConeSet
		global ConeListList
		global HullInfoMap
		if Index == len(ConeListList):
			if len(NewCone.rays()) == 1:
				for Ray in NewCone.rays():
					ConeSet.add(tuple(Ray.list()))
		else:
			for i in xrange(len(ConeListList[Index])):
				TempCone = NewCone.intersection(ConeListList[Index][i])
				if TempCone.dim() > 0:
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
		ConeSetList = []
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

			ConeSetListIndex = len(ConeSetList)
			ConeSetList.append(set())
			PretropEdges = set()
			NotPretropEdges = set()
			ACone = Cone(AEdge.InnerNormals)
			while(len(EdgesToTest) > 0):
				TestEdge = EdgesToTest.pop()
				BEdge = BFaces[0][TestEdge]
				TempCone = Cone(BEdge.InnerNormals).intersection(ACone)
				if len(TempCone.rays()) > 0:
					PretropEdges.add(TestEdge)
					ConeSetList[ConeSetListIndex].add(TempCone)
					for Neighbor in BEdge.Neighbors:
						if Neighbor[0] not in PretropEdges and Neighbor[0] not in NotPretropEdges:
							EdgesToTest.add(Neighbor[0])
				else:
					NotPretropEdges.add(TestEdge)

		global ConeListList
		ConeListList = []
		for i in xrange(len(ConeSetList)):
			ConeListList.append([])
			for NewCone in ConeSetList[i]:
				ConeListList[i].append(NewCone)

		for i in xrange(len(ConeListList[0])):
			IntersectCones(1, ConeListList[0][i])

	ConeList = list(ConeSet)
	ConeList.sort()
	"""
	global Rays
	if len(ConeList) == len(Rays):
		for i in xrange(len(ConeList)):
			if list(ConeList[i]) != Rays[i]:
				print "UNEQUAL!", Rays[i], ConeList[i]
				return 0, 0
	else:
		print "Lists are unequal lengths", len(ConeList), len(Rays)
		return 0, 0
	"""

	print "NewAlg took", time() - NewAlgStart, "seconds."
	print "NewAlg found", len(ConeList), "rays."
	return str(time() - NewAlgStart) + ',' + str(len(ConeList)) # len(ConeList), time() - NewAlgStart#, len(HullInfoMap[(0,"Faces")][0]), len(HullInfoMap[(1,"Faces")][0]), len(HullInfoMap[(2,"Faces")][0])
	
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
	CayleyPts = []

	ToAppend = [[0 for i in xrange(len(PtsList)-1)]]
	for i in xrange(1,len(PtsList)):
		ToAppend.append([1 if i-1 == j else 0 for j in xrange(len(PtsList)-1)])

	for i in xrange(len(PtsList)):
		for j in xrange(len(PtsList[i])):
			CayleyPts.append(PtsList[i][j] + ToAppend[i])
	for i in xrange(len(CayleyPts)):
		for j in xrange(len(CayleyPts[j])):
			CayleyPts[i][j] = Integer(CayleyPts[i][j])

	Faces, IndexToPointMap, PointToIndexMap, Pts = GiftWrap(CayleyPts,True)

	Dim = len(PtsList[0][0])
	PtsMap = {}
	for i in xrange(len(PtsList)):
		for Pt in PtsList[i]:
			PtsMap[tuple(Pt)] = i
	NormalList = []
	for Face in Faces[len(Faces) - 1]:
		CountList = [0 for i in xrange(len(PtsList))]
		for Pt in Face.Vertices:
			CountList[PtsMap[tuple(IndexToPointMap[Pt][0:Dim])]] += 1
		ShouldAddNormal = True
		for i in xrange(len(CountList)):
			if CountList[i] < 2:
				ShouldAddNormal = False
				break
		if ShouldAddNormal == True:
			NormalList.append(Face.InnerNormals[0][0:Dim])

	print "Cayley took", time() - CayleyStart, "seconds."
	print "Cayley found", len(NormalList), "rays."
	return 

#-------------------------------------------------------------------------------
def DoGfan(Polys, R):
	StartTime = time()
	global Rays
	try:
		Rays = R.ideal(Polys).groebner_fan().tropical_intersection().rays()
	except:
		print "Gfan aborted!"
		Rays = []
	for i in xrange(len(Rays)):
		Rays[i] = [-Rays[i][j] for j in xrange(len(Rays[i]))]
	Rays.sort()
	print "Gfan took", time() - StartTime, "seconds."
	print "Gfan found", len(Rays), "rays."
	return

#-------------------------------------------------------------------------------
def DoNaiveAlgorithm(HullInfoMap, NumberOfPolytopes):
	def IntersectCones(Index, NewCone, NumberOfPolytopes):
		global ConeSet
		global HullInfoMap
		Faces = HullInfoMap[(Index + 1,"Faces")]
		for i in xrange(len(Faces[0])):
			TempCone = NewCone.intersection(Cone(Faces[0][i].InnerNormals))
			if TempCone.dim() > 0:
				if Index == NumberOfPolytopes - 2:
					if len(TempCone.rays()) == 1:
						for Ray in TempCone.rays():
							ConeSet.add(tuple(Ray.list()))
				else:
					IntersectCones(Index+1, TempCone, NumberOfPolytopes)
		return

	StartTime = time()
	AFaces = HullInfoMap[(0,"Faces")]
	for i in xrange(len(AFaces[0])):
		IntersectCones(0, Cone(AFaces[0][i].InnerNormals), NumberOfPolytopes)
	
	print "Naive algorithm took", time() - StartTime, "seconds."
	print "Naive algorithm found", len(Rays), "rays."
	return

MyFile = open('testing3','w')
for i in xrange(100):
	MyFile.write(DoTests(3))
	MyFile.write('\n')
MyFile.close()

MyFile = open('testing4','w')
for i in xrange(100):
	MyFile.write(DoTests(4))
	MyFile.write('\n')
MyFile.close()

MyFile = open('testing5','w')
for i in xrange(100):
	MyFile.write(DoTests(5))
	MyFile.write('\n')
MyFile.close()

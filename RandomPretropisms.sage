load("Pretropism_Util.sage")
load("PretropismConfig.sage")
from time import time
import os

#-------------------------------------------------------------------------------
def DoNewAlgorithm(PolysAsPts, HullInfoMap, NumberOfPolytopes):

	def Explore(PolytopeIndex, NewCone):
		BFaces = HullInfoMap[(PolytopeIndex,"Faces")]
		BIndexToPointMap =	HullInfoMap[(PolytopeIndex,"IndexToPointMap")]
		BPointToIndexMap = HullInfoMap[(PolytopeIndex,"PointToIndexMap")]
		B =	HullInfoMap[(PolytopeIndex,"Pts")]

		Normal = [0 for i in xrange(len(NewCone.rays()[0]))]
		for Ray in NewCone.rays():
			Normal = [Normal[i] + Ray[i] for i in xrange(len(Normal))]
		BInitialForm = FindInitialForm(B, Normal)

		ConeList = []
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

		global ConeIntersectionCount
		PretropEdges = set()
		NotPretropEdges = set()
		while(len(EdgesToTest) > 0):
			TestEdge = EdgesToTest.pop()
			BEdge = BFaces[0][TestEdge]
			TempCone = (BEdge.MyCone).intersection(NewCone)
			ConeIntersectionCount += 1
			if len(TempCone.rays()) > 0:
				PretropEdges.add(TestEdge)
				ConeList.append(TempCone)
				for Neighbor in BEdge.Neighbors:
					if Neighbor[0] not in PretropEdges and Neighbor[0] not in NotPretropEdges:
						EdgesToTest.add(Neighbor[0])
			else:
				NotPretropEdges.add(TestEdge)
		return ConeList

	def IntersectCones(Index, NewCone, NumberOfPolytopes):
		global ConeSet
		global HullInfoMap
		if Index == NumberOfPolytopes:
			if len(NewCone.rays()) == 1:
				for Ray in NewCone.rays():
					ConeSet.add(tuple(Ray.list()))
		else:
			CONESS = Explore(Index,NewCone)
			for TempCone in CONESS:
				IntersectCones(Index + 1, TempCone, NumberOfPolytopes)
		return

	NewAlgStart = time()
	global ConeSet
	global ConeIntersectionCount
	ConeIntersectionCount = 0
	ConeSet = set()
	AFaces = HullInfoMap[(0,"Faces")]
	AIndexToPointMap =	HullInfoMap[(0,"IndexToPointMap")]
	APointToIndexMap = HullInfoMap[(0,"PointToIndexMap")]
	A =	HullInfoMap[(0,"Pts")]


	ConeSet = set()
	for AEdge in AFaces[0]:
		IntersectCones(1,Cone(AEdge.InnerNormals), NumberOfPolytopes)

	ConeList = list(ConeSet)
	ConeList.sort()
	return time() - NewAlgStart

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
	return time() - MinkowskiStart

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

	return time() - CayleyStart

#-------------------------------------------------------------------------------
def DoGfan(PolyString, GfanPolys):
	global GfanFileInputPath
	global GfanPath
	f = open(GfanFileInputPath,'w')
	FirstLine = 'Q[ ', PolyString, ' ]'
	f.write('Q[ '+ PolyString+ ' ]')
	f.write('\n')
	SecondLine = '{', str(GfanPolys)[1:-1], '}'
	f.write('{'+ str(GfanPolys)[1:-1]+ '}')
	f.close()
	
	StartGfanTime = time()
	os.system(GfanPath + 'gfan_tropicalintersection < ' + GfanFileInputPath + ' --nocones')
	return time() - StartGfanTime

#-------------------------------------------------------------------------------
def DoGfanFromSage(Polys, R):
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
	return time() - StartTime

#-------------------------------------------------------------------------------
def DoConeIntersectionAlgorithm(HullInfoMap, NumberOfPolytopes):
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
	global ConeSet
	ConeSet = set()
	AFaces = HullInfoMap[(0,"Faces")]
	for i in xrange(len(AFaces[0])):
		IntersectCones(0, Cone(AFaces[0][i].InnerNormals), NumberOfPolytopes)

	return time() - StartTime

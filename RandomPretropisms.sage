load("Pretropism_Util.sage")
load("PretropismConfig.sage")
from time import time
import os

#-------------------------------------------------------------------------------
def DoNewAlgorithm(PolysAsPts, HullInfoMap, NumberOfPolytopes):
	"""
	This is the implementation of our new algorithm to compute pretropisms.
	"""
	def TraverseEdgeSkeleton(PolytopeIndex, NewCone):
		"""
		Explores the edge skeleton and returns a list of cones that are the
		intersection of NewCone with every edge in the pretropism graph.
		"""
		Faces = HullInfoMap[(PolytopeIndex,"Faces")]
		PointToIndexMap = HullInfoMap[(PolytopeIndex,"PointToIndexMap")]
		Pts =	HullInfoMap[(PolytopeIndex,"Pts")]

		Normal = [0 for i in xrange(len(NewCone.rays()[0]))]
		for Ray in NewCone.rays():
			Normal = [Normal[i] + Ray[i] for i in xrange(len(Normal))]
		InitialForm = FindInitialForm(Pts, Normal)

		ConeList = []
		EdgesToTest = set()
		InitialIndices = set([PointToIndexMap[tuple(Pt)] for Pt in InitialForm])
		for i in xrange(len(Faces[0])):
			if len(InitialIndices.intersection(Faces[0][i].Vertices)) != 0:
				EdgesToTest.add(i)

		global ConeIntersectionCount
		PretropGraphEdges = set()
		NotPretropGraphEdges = set()
		while(len(EdgesToTest) > 0):
			TestEdge = EdgesToTest.pop()
			Edge = Faces[0][TestEdge]
			TempCone = (Edge.MyCone).intersection(NewCone)
			ConeIntersectionCount += 1
			if len(TempCone.rays()) > 0:
				PretropGraphEdges.add(TestEdge)
				ConeList.append(TempCone)
				for Neighbor in Edge.Neighbors:
					if Neighbor[0] not in PretropGraphEdges and Neighbor[0] not in NotPretropGraphEdges:
						EdgesToTest.add(Neighbor[0])
			else:
				NotPretropGraphEdges.add(TestEdge)
		return ConeList

	def IntersectCones(Index, NewCone, NumberOfPolytopes):
		"""
		Performs the recursion
		"""
		global ConeSet
		if Index == NumberOfPolytopes:
			if len(NewCone.rays()) == 1:
				for Ray in NewCone.rays():
					ConeSet.add(tuple(Ray.list()))
		else:
			for TempCone in TraverseEdgeSkeleton(Index,NewCone):
				IntersectCones(Index + 1, TempCone, NumberOfPolytopes)
		return

	NewAlgStart = time()
	global ConeSet
	global ConeIntersectionCount
	ConeIntersectionCount = 0
	ConeSet = set()

	for Edge in HullInfoMap[(0,"Faces")][0]:
		IntersectCones(1,Edge.MyCone, NumberOfPolytopes)

	ConeList = list(ConeSet)
	ConeList.sort()
	#print "Number of cone intersections = ", ConeIntersectionCount
	return time() - NewAlgStart

#-------------------------------------------------------------------------------
def DoMinkowskiSum(Polys, PtsList, UseGiftWrap = True):
	"""
	Implements the Minkowski sum method to compute pretropisms.
	"""
	MinkowskiStart = time()
	ProdPoly = 1
	for Poly in Polys:
		ProdPoly = ProdPoly*Poly
	ProdPolysAsPts = [[Integer(j) for j in i] for i in ProdPoly.exponents()]

	if UseGiftWrap:
		Faces, IndexToPointMap, PointToIndexMap, Pts = GiftWrap(ProdPolysAsPts,True)
	else:
		Faces, IndexToPointMap, PointToIndexMap, Pts = WrapHull(ProdPolysAsPts)
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
def DoCayleyPolytope(PtsList, UseGiftWrap = True):
	"""
	Implements the Cayley embedding method to compute pretropisms
	"""

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

	if UseGiftWrap:
		Faces, IndexToPointMap, PointToIndexMap, Pts = GiftWrap(CayleyPts,True)
	else:
		Faces, IndexToPointMap, PointToIndexMap, Pts = WrapHull(CayleyPts)

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
	"""
	Calls the current version of Gfan and computes the tropical intersection. This
	doesn't have the nice wrapping of DoGfanFromSage, but it won't have any wasted
	time making objects. This is only useful for doing direct speed comparisons
	between our code and Gfan.
	"""
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
	"""
	Calls gfan 0.3 via the Sage interface. Sets up the global variable
	Rays which contains the sorted list of pretropisms.
	"""
	StartTime = time()
	global Rays
	try:
		Rays = R.ideal(Polys).groebner_fan().tropical_intersection().rays()
	except:
		print "Gfan aborted!"
		Rays = []
	GfanTime = time() - StartTime
	for i in xrange(len(Rays)):
		Rays[i] = [-Rays[i][j] for j in xrange(len(Rays[i]))]
	Rays.sort()
	return GfanTime

#-------------------------------------------------------------------------------
def DoConeIntersectionAlgorithm(HullInfoMap, NumberOfPolytopes):
	"""
	Function that performs the cone intersection method with pruning along
	the way.
	"""
	def IntersectCones(Index, NewCone, NumberOfPolytopes):
		global ConeSet
		global HullInfoMap
		Faces = HullInfoMap[(Index + 1,"Faces")]
		for i in xrange(len(Faces[0])):
			TempCone = NewCone.intersection(Faces[0][i].MyCone)
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
	for Edge in HullInfoMap[(0,"Faces")][0]:
		IntersectCones(0, Edge.MyCone, NumberOfPolytopes)

	return time() - StartTime

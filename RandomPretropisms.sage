load("Pretropism_Util.sage")
load("PretropismConfig.sage")
from time import time
from multiprocessing import Pool
import os

#-------------------------------------------------------------------------------
def IntersectConesWrapper(ConeIndex):
	global NumberOfPolytopes
	global StartIndex
	ConeIntersectionCount = 0
	ConeContainsCount = 0
	ConeSet = set()
	def TraverseEdgeSkeleton(PolytopeIndex, NewCone, ConeIntersectionCount, ConeContainsCount):
		"""
		Explores the edge skeleton and returns a list of cones that are the
		intersection of NewCone with every edge in the pretropism graph.
		"""
		Faces = HullInfoMap[(PolytopeIndex,"Faces")]
		PointToIndexMap = HullInfoMap[(PolytopeIndex,"PointToIndexMap")]
		Pts = HullInfoMap[(PolytopeIndex,"Pts")]

		Normal = [0 for i in xrange(len(NewCone.rays()[0]))]
		for Ray in NewCone.rays():
			Normal = [Normal[i] + Ray[i] for i in xrange(len(Normal))]
		InitialForm = FindInitialForm(Pts, Normal)

		SkeletonConeSet = set([])
		EdgesToTest = set()
		InitialIndices = set([PointToIndexMap[tuple(Pt)] for Pt in InitialForm])
		for i in xrange(len(Faces[0])):
			if len(InitialIndices.intersection(Faces[0][i].Vertices)) != 0:
				EdgesToTest.add(i)

		PretropGraphEdges = set()
		NotPretropGraphEdges = set()
		while(len(EdgesToTest) > 0):
			TestEdge = EdgesToTest.pop()
			Edge = Faces[0][TestEdge]
			if ConeContains(Edge.MyCone, NewCone):
				ConeContainsCount += 1
				PretropGraphEdges.add(TestEdge)
				SkeletonConeSet.add(NewCone)
				for Neighbor in Edge.Neighbors:
					if Neighbor[0] not in PretropGraphEdges and Neighbor[0] not in NotPretropGraphEdges:
						EdgesToTest.add(Neighbor[0])
			else:
				TempCone = (Edge.MyCone).intersection(NewCone)
				ConeIntersectionCount += 1
				if len(TempCone.rays()) > 0:
					PretropGraphEdges.add(TestEdge)
					SkeletonConeSet.add(TempCone)
					for Neighbor in Edge.Neighbors:
						if Neighbor[0] not in PretropGraphEdges and Neighbor[0] not in NotPretropGraphEdges:
							EdgesToTest.add(Neighbor[0])
				else:
					NotPretropGraphEdges.add(TestEdge)
		return SkeletonConeSet, ConeIntersectionCount, ConeContainsCount

	def IntersectCones(Index, NewCone):
		"""
		Performs the recursion
		"""
		if (Index != NumberOfPolytopes) and (Index == StartIndex + 1):
			if len(NewCone.rays()) > 0:
				ConeSet.add(NewCone)
			return 0, 0
		elif Index == NumberOfPolytopes:
			if len(NewCone.rays()) > 0:
				for Ray in NewCone.rays():
					ConeSet.add(tuple(Ray.list()))
			return 0, 0
		else:
			SkeletonConeSet, NewConeIntersectionCount, NewConeContainsCount = TraverseEdgeSkeleton(Index, NewCone, 0, 0)
			ConeIntersectionCount = NewConeIntersectionCount
			ConeContainsCount = NewConeContainsCount
			for TempCone in SkeletonConeSet:
				NewConeIntersectionCount, NewConeContainsCount = IntersectCones(Index + 1, TempCone)
				ConeIntersectionCount += NewConeIntersectionCount
				ConeContainsCount += NewConeContainsCount
		return ConeIntersectionCount, ConeContainsCount
	MyStart = time()
	global ConesToFollowList
	return ConeSet, IntersectCones(StartIndex, ConesToFollowList[ConeIndex]), time() - MyStart

#-------------------------------------------------------------------------------
def DoNewAlgorithm(HullInfoMaps, ThreadCount):
	"""
	This is the implementation of our new algorithm to compute pretropisms.
	"""
	
	ConeIntersectionCount = 0
	ConeContainsCount = 0
	NewAlgStart = time()
	ConesToFollow = set()
	global StartIndex
	global NumberOfPolytopes
	global ConesToFollowList
	ConesToFollowList = [Edge.MyCone for Edge in HullInfoMap[(0,"Faces")][0]]
	for i in xrange(1,NumberOfPolytopes):
		print ""
		print i, "out of", NumberOfPolytopes-1
		StartIndex = i
		MyPool = Pool(ThreadCount)
		ResultList = MyPool.map(IntersectConesWrapper, [i for i in xrange(len(ConesToFollowList))])
		#print ResultList
		MyPool.terminate()
		
		ConesToFollowSet = set()
		TimeList = []
		for Element in ResultList:
			ConesToFollowSet = ConesToFollowSet.union(Element[0])
			ConeIntersectionCount += Element[1][0]
			ConeContainsCount += Element[1][1]
			TimeList.append(Element[2])
		TimeList.sort()
		#print "List of times", TimeList
		ConesToFollowList = list(ConesToFollowSet)
		print len(ConesToFollowList), len(ConesToFollowSet), ConeIntersectionCount, ConeContainsCount, TimeList[0:5], TimeList[len(TimeList)-6:len(TimeList)-1]

	ConesToFollowList.sort()
	TimeList.sort()
	#print "List of times", TimeList
	print "NEW RESULT", ConesToFollowList
	print "Number of cone intersections = ", ConeIntersectionCount
	print "Number of cone contains = ", ConeContainsCount
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
		Rays[i] = tuple(-Rays[i][j] for j in xrange(len(Rays[i])))
	Rays.sort()
	print "GFANRESULT", Rays
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

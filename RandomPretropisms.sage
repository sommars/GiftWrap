load("Pretropism_Util.sage")
load("PretropismConfig.sage")
from time import time
from multiprocessing import Pool
import os
from sage.libs.ppl import Variable, Constraint_System, C_Polyhedron
#-------------------------------------------------------------------------------
def IntersectConesWrapper(ConeIndex):

	global NumberOfPolytopes
	global StartIndex
	global InFinalRound
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
	Counts = list(IntersectCones(StartIndex, ConesToFollowList[ConeIndex]))

	if StartIndex != NumberOfPolytopes - 1:
		NewTuples = []
		for MyCone in list(ConeSet):
			NewTuples.append((MyCone, GetCPolyhedron(MyCone)))
		NewCones = set([])
		ContainCount = 0
		IndicesToIgnore = set()
		for j in xrange(len(NewTuples)):
			TempCone = NewTuples[j]
			ConsSystem = TempCone[1]
			ShouldAddCone = True
			MyDim = TempCone[0].dim()
			for k in xrange(len(NewTuples)):
				if j == k or k in IndicesToIgnore:
					continue
				TestCone = NewTuples[k]
				if TestCone[0].dim() < MyDim:
					continue
				ContainCount += 1
				if TestCone[1].contains(ConsSystem):
					IndicesToIgnore.add(j)
					ShouldAddCone = False
					break
			if ShouldAddCone == True:
				NewCones.add(TempCone[0])
		Counts[1] = Counts[1] + ContainCount
	else:
		NewCones = set([])
		for MyCone in list(ConeSet):
			NewCones.add(MyCone)
	return NewCones, Counts, time() - MyStart

#-------------------------------------------------------------------------------
def DoNewAlgorithm(HullInfoMaps, ThreadCount):
	"""
	This is the implementation of our new algorithm to compute pretropisms.
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
		ResultList = MyPool.map(IntersectConesWrapper, [j for j in xrange(len(ConesToFollowList))])
		MyPool.terminate()
		
		ConesToFollowSet = set()
		TimeList = []
		for Element in ResultList:
			ConesToFollowSet = ConesToFollowSet.union(Element[0])
			ConeIntersectionCount += Element[1][0]
			ConeContainsCount += Element[1][1]
			TimeList.append(Element[2])
		TimeList.sort()
		ConesToFollowList = list(ConesToFollowSet)
		print len(ConesToFollowList), len(ConesToFollowSet), ConeIntersectionCount, ConeContainsCount, TimeList[0:5], TimeList[len(TimeList)-6:len(TimeList)-1]
		NewCones = set()
		if i == NumberOfPolytopes - 1:
			break

		if True == True:
			NewTuples = []
			IndicesToIgnore = set()
			for MyCone in ConesToFollowList:
				NewTuples.append((MyCone, GetCPolyhedron(MyCone)))
			MyStartTime = time()
			for j in xrange(len(NewTuples)):
				TempCone = NewTuples[j]
				ConsSystem = TempCone[1]
				ShouldAddCone = True
				MyDim = TempCone[0].dim()
				for k in xrange(len(NewTuples)):
					if j == k or k in IndicesToIgnore:
						continue
					TestCone = NewTuples[k]
					if TestCone[0].dim() < MyDim:
						continue
					ConeContainsCount += 1
					if TestCone[1].contains(ConsSystem):
						IndicesToIgnore.add(j)
						ShouldAddCone = False
						break
				if ShouldAddCone == True:
					NewCones.add(TempCone[0])

			print "JEFF", len(NewCones), len(ConesToFollowList), time() - MyStartTime
			ConesToFollowList = list(NewCones)

	ConesToFollowList.sort()
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

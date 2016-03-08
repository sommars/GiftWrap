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
	ConeIntersectionCount = 0
	ConeContainsCount = 0
	def TraverseEdgeSkeleton(PolytopeIndex, InputReducedCone):
		"""
		Explores the edge skeleton and returns a list of cones that are the
		intersection of NewCone with every edge in the pretropism graph.
		"""
		ConeIntersectionCount = 0
		ConeContainsCount = 0
		Faces = HullInfoMap[(PolytopeIndex,"Faces")]
		PointToIndexMap = HullInfoMap[(PolytopeIndex,"PointToIndexMap")]
		Pts = HullInfoMap[(PolytopeIndex,"Pts")]

		Normal = [0 for i in xrange(len(InputReducedCone.MyCone.rays()[0]))]
		for Ray in InputReducedCone.MyCone.rays():
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
		if not InputReducedCone.Has_CPolyhedron:
			InputReducedCone.CPolyhedron = GetCPolyhedron(InputReducedCone.MyCone)
			InputReducedCone.Has_CPolyhedron = True

		while(len(EdgesToTest) > 0):
			TestEdgeIndex = EdgesToTest.pop()
			Edge = Faces[0][TestEdgeIndex]
			if Edge.CPolyhedron.contains(InputReducedCone.CPolyhedron):
				ConeContainsCount += 1
				PretropGraphEdges.add(TestEdgeIndex)
				SkeletonConeSet.add(InputReducedCone)
				for Neighbor in Edge.Neighbors:
					if Neighbor[0] not in PretropGraphEdges and Neighbor[0] not in NotPretropGraphEdges:
						EdgesToTest.add(Neighbor[0])
			else:
				#Here's where the dimension reduction stuff will need to happen.
				global TempEquations
				TempEquations = deepcopy(Edge.Equations)
				global TempInequalities
				TempInequalities = deepcopy(Edge.Inequalities)
				for Reduction in InputReducedCone.Reductions:
					for j in xrange(len(TempEquations)):
						TempEquations[j] = DropDim(Reduction[1],TempEquations[j],Reduction[0])
					for j in xrange(len(TempInequalities)):
						TempInequalities[j] = DropDim(Reduction[1],TempInequalities[j],Reduction[0])
				Poly = Polyhedron(eqns=TempEquations,ieqs=TempInequalities)
				if len(Poly.rays_list()) + len(Poly.lines_list()) == 0:
					NotPretropGraphEdges.add(TestEdgeIndex)
					continue
				TempReducedCone = Cone(Poly.rays_list() + Poly.lines_list() + [[-Coord for Coord in TempLine ] for TempLine in Poly.lines_list()])
				TempCone = TempReducedCone.intersection(InputReducedCone.ReducedCone)
				ConeIntersectionCount += 1
				if len(TempCone.rays()) > 0:
					PretropGraphEdges.add(TestEdgeIndex)
					NewReducedCone = ReducedConeClass()
					NewReducedCone.EdgeTuples = InputReducedCone.EdgeTuples + [(PolytopeIndex, TestEdgeIndex)]
					NewReducedCone.MyCone = (Edge.MyCone).intersection(InputReducedCone.MyCone)
					if NewReducedCone.MyCone.dim() != TempCone.dim():
						print "This is bad"
						print avariablethatdoesntexist
					NewReducedCone.Reductions = deepcopy(InputReducedCone.Reductions)
					#Now we need to check if the new cone is reducible
					if TempCone.lattice_dim() == TempCone.dim():
						NewReducedCone.ReducedCone = TempCone
					else:
						TempPolyhedron = TempCone.polyhedron()
						TempEquations = TempPolyhedron.equations_list()
						TempInequalities = TempPolyhedron.inequalities_list()
						for j in xrange(len(TempEquations)):
							Equation = TempEquations[j]
							Index = next((i for i, x in enumerate(Equation) if x), None)
							NewReducedCone.Reductions.append([Index, TempEquations[j]])
							for k in xrange(j+1, len(TempEquations)):
								TempEquations[k] = DropDim(Equation,TempEquations[k],Index)
							for k in xrange(len(TempInequalities)):
								TempInequalities[k] = DropDim(Equation,TempInequalities[k],Index)
						TempPolyhedron = Polyhedron(ieqs=TempInequalities)
						NewReducedCone.ReducedCone = Cone(TempPolyhedron.rays_list() + TempPolyhedron.lines_list() + [[-Coord for Coord in TempLine ] for TempLine in TempPolyhedron.lines_list()])
					SkeletonConeSet.add(NewReducedCone)
					for Neighbor in Edge.Neighbors:
						if Neighbor[0] not in PretropGraphEdges and Neighbor[0] not in NotPretropGraphEdges:
							EdgesToTest.add(Neighbor[0])
				else:
					NotPretropGraphEdges.add(TestEdgeIndex)
		return SkeletonConeSet, ConeIntersectionCount, ConeContainsCount

	def GetCPolyhedron(OldCone):
		cs = Constraint_System()
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

	def DropDim(A,B,Index):
		if B[Index] == 0:
			del B[Index]
			return B
		AMult = lcm(A[Index],B[Index])/abs(A[Index])*sign(B[Index])
		BMult = lcm(A[Index],B[Index])/abs(B[Index])*sign(A[Index])
		NewB = []
		for i in xrange(len(B)):
			if i == Index:
				continue
			NewB.append(-AMult*A[i]+BMult*B[i])
		return NewB

	MyStart = time()
	global ConesToFollowList
	InputReducedCone = ConesToFollowList[ConeIndex]
	SkeletonReducedConeSet, ConeIntersectionCount, ConeContainsCount = TraverseEdgeSkeleton(StartIndex, ConesToFollowList[ConeIndex])
	Counts = [ConeIntersectionCount, ConeContainsCount]

	ConeList = list(SkeletonReducedConeSet)
	for ReducedCone in ConeList:
		if not ReducedCone.Has_CPolyhedron:
			ReducedCone.CPolyhedron = GetCPolyhedron(ReducedCone.MyCone)
			ReducedCone.Has_CPolyhedron = True
	NewCones = set([])
	ContainCount = 0
	IndicesToIgnore = set()
	for j in xrange(len(ConeList)):
		ReducedCone = ConeList[j]
		ShouldAddCone = True
		MyDim = ReducedCone.MyCone.dim()
		for k in xrange(len(ConeList)):
			ReducedCone2 = ConeList[k]
			if j == k or k in IndicesToIgnore:
				continue
			if ReducedCone2.MyCone.dim() < MyDim:
				continue
			ContainCount += 1
			if ReducedCone2.CPolyhedron.contains(ReducedCone.CPolyhedron):
				IndicesToIgnore.add(j)
				ShouldAddCone = False
				break

		if ShouldAddCone == True:
			NewCones.add(ReducedCone)
	Counts[1] = Counts[1] + ContainCount
	return NewCones, Counts, time() - MyStart

#-------------------------------------------------------------------------------
def DoConeContainment(Index):
	global ConesToFollowList
	global ConeConstraintSets
	global ConstraintList
	ConstraintTestMap = {}
	ReducedCone = ConesToFollowList[Index]
	TestCPolyhedron = ReducedCone.CPolyhedron
	ShouldAddCone = True
	MyDim = ReducedCone.MyCone.dim()
	for k in xrange(len(ConesToFollowList)):
		if Index == k:
			continue

		TestCone = ConesToFollowList[k]
		if TestCone.MyCone.dim() < MyDim:
			continue
		
		TestConstraintList = list(ConeConstraintSets[k])
		for ii in xrange(len(TestConstraintList)):
			BoolValue = TestConstraintList[ii]
			if not ConstraintTestMap.has_key(BoolValue):
				if ConstraintList[BoolValue].contains(TestCPolyhedron):
					ConstraintTestMap[BoolValue] = 1
				else:
					ConstraintTestMap[BoolValue] = 0
			if ConstraintTestMap[BoolValue] != 1:
				break
			if ii == len(TestConstraintList) - 1:
				ShouldAddCone = False
				return -1
	return ReducedCone

#-------------------------------------------------------------------------------
def DoNewAlgorithm(HullInfoMaps, ThreadCount):
	"""
	This is the implementation of our new algorithm to compute pretropisms.
	"""

	def DropDim(A,B,Index):
		if B[Index] == 0:
			del B[Index]
			return B
		AMult = lcm(A[Index],B[Index])/abs(A[Index])*sign(B[Index])
		BMult = lcm(A[Index],B[Index])/abs(B[Index])*sign(A[Index])
		NewB = []
		for i in xrange(len(B)):
			if i == Index:
				continue
			NewB.append(-AMult*A[i]+BMult*B[i])
		return NewB

	ConeIntersectionCount = 0
	ConeContainsCount = 0
	NewAlgStart = time()
	ConesToFollow = set()
	global StartIndex
	global NumberOfPolytopes
	global ConesToFollowList
	ConesToFollowList = []
	for i in xrange(len(HullInfoMap[(0,"Faces")][0])):
		NewFace = HullInfoMap[(0,"Faces")][0][i]
		NewReducedCone = ReducedConeClass()
		NewReducedCone.EdgeTuples.append((0,i))
		NewReducedCone.MyCone = NewFace.MyCone
		TempEquations = deepcopy(NewFace.Equations)
		TempInequalities = deepcopy(NewFace.Inequalities)
		for j in xrange(len(TempEquations)):
			Equation = TempEquations[j]
			Index = next((i for i, x in enumerate(Equation) if x), None)
			#add the new equation and the relevant index to the Reductions
			NewReducedCone.Reductions.append([Index, TempEquations[j]])
			#reduce all of the equations in TempEquations
			for k in xrange(j+1, len(TempEquations)):
				TempEquations[k] = DropDim(Equation,TempEquations[k],Index)
			for k in xrange(len(TempInequalities)):
				TempInequalities[k] = DropDim(Equation,TempInequalities[k],Index)

		TempPolyhedron = Polyhedron(ieqs=TempInequalities)
		NewReducedCone.ReducedCone = Cone(TempPolyhedron.rays_list() + TempPolyhedron.lines_list() + [[-Coord for Coord in TempLine ] for TempLine in TempPolyhedron.lines_list()])
		print NewReducedCone.Reductions
		L = []
		for TempRay in NewReducedCone.ReducedCone.rays():
			L.append(list(TempRay))
		LL = []
		for TempRay in NewReducedCone.MyCone.rays():
			LL.append(list(TempRay))
		for TempReduction in NewReducedCone.Reductions:
			for TempRay in LL:
				del TempRay[TempReduction[0]-1]
		L.sort()
		LL.sort()
		ConesToFollowList.append(NewReducedCone)
		if L != LL:
			print "This is bad"
			print avariablethatdoesntexist

	for i in xrange(1, NumberOfPolytopes):
		print ""
		print i, "out of", NumberOfPolytopes-1
		StartIndex = i
		ParallelStartTime = time()
		MyPool = Pool(ThreadCount)
		ResultList = MyPool.map(IntersectConesWrapper, [j for j in xrange(len(ConesToFollowList))])
		MyPool.terminate()
		print "ParallelRealTime", time() - ParallelStartTime
		ParallelCleanUpTime = time()
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
		print "ParallelCleanupTime", time() - ParallelCleanUpTime
		NewCones = set()
		if i == NumberOfPolytopes - 1:
			break

		ConstraintMap = {}
		global ConstraintList
		ConstraintList = []
		ConstraintCount = 0
		global ConeConstraintSets
		ConeConstraintSets = []
		NewTime = time()
		CPolyTime = 0
		ConstrTime = 0
		for ReducedCone in ConesToFollowList:
			TestTime = time()
			NewCPolyhedron = ReducedCone.CPolyhedron
			CPolyTime += time() - TestTime
			NewSet = set()
			ConstrTimeTime = time()
			for Cons in NewCPolyhedron.constraints():
				ConsStr = str(Cons)
				if not ConstraintMap.has_key(ConsStr):
					ConstraintMap[ConsStr] = ConstraintCount
					ConstraintList.append(C_Polyhedron(Constraint_System(Cons)))
					ConstraintCount += 1
				NewSet.add(ConstraintMap[ConsStr])
			ConeConstraintSets.append(NewSet)
			ConstrTime += time() - ConstrTimeTime
		print "ConstraintTime", ConstrTime
		print "CPolyTime", CPolyTime
		print "Number of constraints", len(ConstraintMap)
		print "Upper bound on number of tests", len(ConstraintMap)*len(ConesToFollowList)
		print "Cone contain setup time", time() - NewTime
		MyStartTime = time()
		MyPool = Pool(ThreadCount)
		ResultList = MyPool.map(DoConeContainment, [j for j in xrange(len(ConesToFollowList))])
		MyPool.terminate()
		for NewCone in ResultList:
			if NewCone != -1:
				NewCones.add(NewCone)

		print "CONTAINMENT END", len(NewCones), len(ConesToFollowList), time() - MyStartTime
		ConesToFollowList = list(NewCones)
		
	print "GENERATING RESULTS..."
	Result = set([])
	for ReducedCone in ConesToFollowList:
		EdgeTuples = ReducedCone.EdgeTuples
		if len(EdgeTuples) == 0:
			continue
		NewCone = HullInfoMap[(EdgeTuples[0][0],"Faces")][0][EdgeTuples[0][1]].MyCone
		for i in xrange(1, len(EdgeTuples)):
			NewCone = NewCone.intersection(HullInfoMap[(EdgeTuples[i][0],"Faces")][0][EdgeTuples[i][1]].MyCone)
		for Ray in NewCone.rays():
			Result.add(tuple(Ray.list()))
			
	Result = list(Result)
	Result.sort()
	print "NEW RESULT", Result
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

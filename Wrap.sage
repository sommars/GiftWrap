from time import time
load("GiftUtil.sage")
# The code in this unit is the code that will have to change for going between 3 and n dimensions.

#If all of the points have a column that's identical, this will break. I should immediately make a map in that case..?
def GiftWrap(Pts):
	global Faces
	global InitialDim
	global PointToIndexMap
	global IndexToPointMap
	if not PtsAreValid(Pts):
		print "The input set of points is not valid."
		raw_input()
		return
	Pts = RemoveDups(Pts)
	OriginalBarycenter = FindBarycenter(Pts)
	if PointToIndexMap == "Start":
		PointToIndexMap, IndexToPointMap = MakeIndexMaps(Pts)
	Pts, ShortPointToLongPointMap, LongPointToShortPointMap, LocalDim = WrapMaps(Pts)
	if Faces == "Start":
		Faces = []
		InitialDim = LocalDim
		for Dim in xrange(1,len(Pts[0])):
			Faces.append([])
	Barycenter = FindBarycenter(Pts)
	if InitialDim <= 2:
		print "Still need to handle low dimensional situations gracefully"
		return 
	FirstPts = FindInitialFacet(Pts, Barycenter)
	IndexOfFirstFace = MakeFace(FirstPts, Pts, LongPointToShortPointMap, ShortPointToLongPointMap, Barycenter)
	FaceIndices = [IndexOfFirstFace]
	Counter = 0
	EdgesToRotateOver = set([])
	for ChildEdge in Faces[LocalDim-2][IndexOfFirstFace].Children:
		EdgesToRotateOver.add((IndexOfFirstFace, ChildEdge))
	while len(EdgesToRotateOver) > 0:
		Edge = EdgesToRotateOver.pop()
		FaceWithEdge = Faces[LocalDim-2][Edge[0]]
		ChildFace = Faces[LocalDim-3][Edge[1]]

		Counter += 1
		if Counter == 200:
			print "Internal error, computation time exceeded"
			print len(Faces[0])
			print len(Faces[1])
			raw_input()

		VerticesToUse = []
		for Vertex in FaceWithEdge.Vertices:
			VerticesToUse.append(LongPointToShortPointMap[tuple(IndexToPointMap[Vertex])])

		EdgePts = []
		for Vertex in ChildFace.Vertices:
			EdgePts.append(LongPointToShortPointMap[tuple(IndexToPointMap[Vertex])])
		Normal = GetNormalFromHNF(GetHNF(VerticesToUse))

		# Need to make the normal an inner normal
		if not NormalPointsTowardsPt(Normal, Barycenter, EdgePts[0]):
			Normal = [-Normal[i] for i in xrange(len(Normal))]
		InnerNormal = Normal
		FacePts = EdgePts + FindNewFacetPtsFromEdge(Pts, EdgePts, InnerNormal, VerticesToUse)[0]
		#Call MakeFacet on the points I've found and do some bookkeeping
		NewFaceIndex = MakeFace(FacePts, Pts, LongPointToShortPointMap, ShortPointToLongPointMap, Barycenter)
		if NewFaceIndex not in FaceIndices:
			EdgesToRotateOver.add(Edge)
			FaceIndices.append(NewFaceIndex)

			#Note that a Neighbor = (Neighbor, ChildFaceShared)
			for ChildEdge in Faces[LocalDim-2][NewFaceIndex].Children:
				ShouldBreak = False
				for Index in FaceIndices:
					if ShouldBreak == True:
						break
					Face = Faces[LocalDim-2][Index]
					for Neighbor in Face.Neighbors:
						if (Neighbor[1] == ChildEdge) and (Neighbor[0] == NewFaceIndex):
							if (Index, ChildEdge) in EdgesToRotateOver:
								EdgesToRotateOver.remove((Index, ChildEdge))
								ShouldBreak = True
								break
				if ShouldBreak == False:
					EdgesToRotateOver.add((NewFaceIndex,ChildEdge))
	"""
	if LocalDim == 3:
		MyFacetPtRefs = set([])
		for Index in FaceIndices:
			MyFacetPtRefs = MyFacetPtRefs.union(Faces[LocalDim-2][Index].Vertices)
		MyFacetPts = []
		for Pt in MyFacetPtRefs:
			MyFacetPts.append(IndexToPointMap[Pt])
		AllPts = [[47, 13, 48, 92], [8, 25, 69, 8], [85, 25, 75, 5], [86, 67, 96, 64], [95, 58, 89, 73], [17, 0, 35, 31]]
		CheckFacetAgainstSage(AllPts, MyFacetPts)
	"""
	Vertices = set([])
	for Index in FaceIndices:
		Vertices = Vertices.union(Faces[LocalDim - 2][Index].Vertices)
	if LocalDim%2 == 0:
		Goal = 0
	else:
		Goal = 2

	FaceSum = 0
	SubFaces = set(FaceIndices)
	LenList = [0 for i in xrange(LocalDim - 1)]
	for MyDim in xrange(LocalDim - 1):
		NextSubFaces = set([])
		CurDim = LocalDim - MyDim - 2
		for Index in SubFaces:
			NextSubFaces = NextSubFaces.union(Faces[CurDim][Index].Children)
		FaceSum += (-1)^(CurDim+1)*len(SubFaces)
		LenList[CurDim] = (-1)^(CurDim+1)*len(SubFaces)
		SubFaces = NextSubFaces
	if len(Vertices) + FaceSum != Goal:
		print "Euler characteristic failed for polytope that lives in dimension", LocalDim
		print "len(Vertices) = ", len(Vertices)
		for i in xrange(len(Faces)):
			print "Length of Faces[", i, "] is equal to", LenList[i]
		print "Goal", Goal
		for i in xrange(len(Faces)):
			print len(Faces[i])
		raw_input()
	if InitialDim == LocalDim:
		#PrintFaces()
		print "Euler characteristic passed for polytope that lives in dimension", LocalDim
		return 
	return FaceIndices

#-------------------------------------------------------------------------------
def Make2dFace(FacePts, Pts, LongPointToShortPointMap, ShortPointToLongPointMap, Barycenter):
	FacePts, LocalShortPointToLongPointMap, LocalLongPointToShortPointMap, LocalDim = WrapMaps(FacePts)
	ShortHull = ConvexHull2d(FacePts)
	Hull = []
	for Pt in ShortHull:
		Hull.append(LocalShortPointToLongPointMap[tuple(Pt)])

	NewFace = Face()
	Normal = GetNormalFromHNF(GetHNF(Hull))
	# Need to make the normal an inner normal
	if not NormalPointsTowardsPt(Normal, Barycenter, Hull[0]):
		Normal = [-Normal[i] for i in xrange(len(Normal))]
	NewFace.InnerNormal = Normal
	Edges = [[Hull[len(Hull) - 1], Hull[0]]]
	for i in xrange(len(Hull) - 1):
		Edges.append([Hull[i],Hull[i+1]])

	NewFace.Vertices = set()
	for Pt in Hull:
		NewFace.Vertices.add(PointToIndexMap[tuple(Pt)])

	Children = set([])
	for i in xrange(len(Hull)):
		EdgeFound = False
		Pt1 = PointToIndexMap[tuple(Hull[i])]
		Pt2 = PointToIndexMap[tuple(Hull[(i+1)%len(Hull)])]
		Edge = set([Pt1, Pt2])
		for i in xrange(len(Faces[0])):
			if Faces[0][i].Vertices == Edge:
				Children.add(i)
				EdgeFound = True
				break
		if EdgeFound == False:
			EdgeFace = Face()
			EdgeFace.Dimension = 1
			EdgeFace.Vertices = Edge
			Faces[0].append(EdgeFace)
			Children.add(len(Faces[0])-1)
			EdgeFace.Neighbors = set([])
			for i in xrange(len(Faces[0])):
				PossibleNeighbor = Faces[0][i]
				if len(PossibleNeighbor.Vertices.intersection(EdgeFace.Vertices)) == 1:
					PossibleNeighbor.Neighbors.add((len(Faces[0])-1,list(PossibleNeighbor.Vertices.intersection(EdgeFace.Vertices))[0]))
					EdgeFace.Neighbors.add((i,list(PossibleNeighbor.Vertices.intersection(EdgeFace.Vertices))[0]))

	NewFace.Children = Children
	NewFace.Dimension = 2
	NewFace.Neighbors = set([])
	LocationOfNewFace = len(Faces[1])
	for i in xrange(len(Faces[1])):	
		PossibleNeighbor = Faces[1][i]
		Intersection = list(PossibleNeighbor.Children.intersection(NewFace.Children))
		if len(Intersection) == 1:
			PossibleNeighbor.Neighbors.add((LocationOfNewFace,Intersection[0]))
			NewFace.Neighbors.add((i,Intersection[0]))
		if len(Intersection) > 1:
			raw_input()
	Faces[1].append(NewFace)
	return LocationOfNewFace

#-------------------------------------------------------------------------------
def PrintFaces():
	print "BEGIN"
	for i in xrange(len(Faces)):
		print ""
		print ""
		print ""
		print "----------------------------------------------"
		print "Below are faces of dim", i + 1
		for Face in Faces[i]:
			Face.PrintProps()
	print "END"
	return

#-------------------------------------------------------------------------------
def PrintFaceLens():
	print "BEGIN"
	for i in xrange(len(Faces)):
		print "LenFaces[",i,"] = ", len(Faces[i])
	print "END"
	return

#-------------------------------------------------------------------------------
def MakeFace(FacePts, Pts, LongPointToShortPointMap, ShortPointToLongPointMap, Barycenter):
	Dimension = CheckNormalFormDim(GetHNF(FacePts))

	#This is to check if it's already been found
	FacePtsInCorrectDim = []
	PtIndices = []
	for Pt in FacePts:
		LongPt = ShortPointToLongPointMap[tuple(Pt)]
		FacePtsInCorrectDim.append(LongPt)
		PtIndices.append(PointToIndexMap[tuple(LongPt)])
	for i in xrange(len(Faces[Dimension-1])):
		if len((Faces[Dimension-1][i].Vertices).difference(set(PtIndices))) == 0:
			return i

	if Dimension == 2:
		return Make2dFace(FacePtsInCorrectDim, Pts, LongPointToShortPointMap, ShortPointToLongPointMap, Barycenter)
	elif Dimension > 2:
		NewFace = Face()
		NewFaceIndex = len(Faces[Dimension-1])
		NewFace.Children = set(GiftWrap(FacePtsInCorrectDim))
		NewFace.Vertices = set([])
		for ChildFaceIndex in NewFace.Children:
			ChildFace = Faces[Dimension-2][ChildFaceIndex]
			NewFace.Vertices = NewFace.Vertices.union(ChildFace.Vertices)
			ChildFace.Parents.append(NewFaceIndex)


		VerticesToUse = []
		for Vertex in NewFace.Vertices:
			VerticesToUse.append(IndexToPointMap[Vertex])

		Normal = GetNormalFromHNF(GetHNF(VerticesToUse))
		# Need to make the normal an inner normal
		if not NormalPointsTowardsPt(Normal, Barycenter, VerticesToUse[0]):
			Normal = [-Normal[i] for i in xrange(len(Normal))]
		NewFace.InnerNormal = Normal
		NewFace.Dimension = Dimension

		# Not entirely sure about this. Simplex vs. Nonsimplex case..?
		NewFace.Neighbors = set([])
		for i in xrange(len(Faces[Dimension - 1])):
			PossibleNeighbor = Faces[Dimension - 1][i]
			Intersection = list(PossibleNeighbor.Children.intersection(NewFace.Children))
			if len(Intersection) == 1:
				PossibleNeighbor.Neighbors.add((NewFaceIndex,Intersection[0]))
				NewFace.Neighbors.add((i,Intersection[0]))
			if len(Intersection) > 1:
				raw_input()
		Faces[Dimension-1].append(NewFace)
		IndexOfFace = NewFaceIndex
	else:
		print "Internal error, dimension ", Dimension, " not expected in MakeFace."
		raw_input()
	return IndexOfFace


def MakeRandomPointSet(Dim,Num):
	Pts = []
	for i in xrange(Num):
		Pts.append([])
		for j in xrange(Dim):
			Pts[i].append(Integer(randint(0,100)))
	return Pts

Tests = []
for i in xrange(10):
	# First param is dimension, second is number of points
	Tests.append((MakeRandomPointSet(5,10), 'Random'))


Cyclic = CreateCyclicLists(4)
Tests = []
for i in xrange(len(Cyclic)):
	Tests.append((Cyclic[i], 'Cyclic'))


for k in range(20):
	print ""
for Test in Tests:
	global Faces
	Faces = "Start"
	global PointToIndexMap
	PointToIndexMap = "Start"
	global InitialDim
	InitialDim = "Sentinel"
	Pts = Test[0]
	print Test[1]
	print Pts
	StartTime = time()
	GiftWrap(Pts)
	EndTime = time()
	print "Done with " + Test[1]
	print "Time = ", EndTime - StartTime
	print ""
	print ""

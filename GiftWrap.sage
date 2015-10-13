from time import time
load("GiftWrap_Util.sage")

def GiftWrap(Pts):
	global Faces
	global PointToIndexMap
	global IndexToPointMap
	global InitialDim
	if not PtsAreValid(Pts):
		print "The input set of points is not valid."
		raw_input()
		return
	Pts = RemoveDups(Pts)
	if len(Pts) < 2:
		print "The point on the convex hull:"
		print Pts[0]
		return
	if PointToIndexMap == "Start":
		PointToIndexMap, IndexToPointMap = MakeIndexMaps(Pts)
	Pts, ShortPointToLongPointMap, LongPointToShortPointMap, LocalDim = WrapMaps(Pts)
	if Faces == "Start":
		Faces = []
		InitialDim = LocalDim
		for Dim in xrange(1,len(Pts[0])):
			Faces.append([])
	Barycenter = FindBarycenter(Pts)
	if InitialDim == 2:
		print "The points on the convex hull in the order in which they occur:"
		FullDimPts = []
		for Pt in ConvexHull2d(Pts):
			FullDimPts.append(ShortPointToLongPointMap[tuple(Pt)])
		print FullDimPts
		return
	elif InitialDim == 1:
		print "The points on the convex hull in the order in which they occur:"
		FullDimPts = [ShortPointToLongPointMap[tuple(Pts[i])] for i in xrange(len(Pts))]
		FullDimPts.sort()
		print FullDimPts[0], FullDimPts[1]
		return
	FirstPts = FindInitialFacet(Pts)
	IndexOfFirstFace = MakeFace(FirstPts, Pts, LongPointToShortPointMap, ShortPointToLongPointMap, Barycenter)
	FaceIndices = set([IndexOfFirstFace])
	EdgesToRotateOver = set([])
	for ChildEdge in Faces[LocalDim-2][IndexOfFirstFace].Children:
		EdgesToRotateOver.add((IndexOfFirstFace, ChildEdge))
	Counter = 0
	MaxCounter = 400
	#Something like this: MaxCounter = len(Combinations(len(Pts),LocalDim).list())
	while len(EdgesToRotateOver) > 0:
		if Counter == MaxCounter:
			print "Internal error, computation time exceeded"
			raw_input()
		
		Edge = EdgesToRotateOver.pop()
		FaceWithEdge = Faces[LocalDim-2][Edge[0]]
		ChildFace = Faces[LocalDim-3][Edge[1]]

		VerticesToUse = []
		for Vertex in FaceWithEdge.Vertices:
			VerticesToUse.append(LongPointToShortPointMap[tuple(IndexToPointMap[Vertex])])

		EdgePts = []
		for Vertex in ChildFace.Vertices:
			EdgePts.append(LongPointToShortPointMap[tuple(IndexToPointMap[Vertex])])

		Normal = GetNormalFromHNF(GetHNF(VerticesToUse))
		InnerNormal = MakeNormalPointInDirectionOfPt(Normal, Barycenter, EdgePts[0])
		FacePts = EdgePts + FindNewFacePtsFromEdge(Pts, EdgePts, InnerNormal, VerticesToUse)[0]

		#Call MakeFace on the points I've found and do some bookkeeping
		NewFaceIndex = MakeFace(FacePts, Pts, LongPointToShortPointMap, ShortPointToLongPointMap, Barycenter)
		if NewFaceIndex not in FaceIndices:
			EdgesToRotateOver.add(Edge)
			FaceIndices.add(NewFaceIndex)

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
		Counter += 1

	# Compute the Euler characteristic to check for correctness
	if InitialDim == LocalDim:
		Vertices = set([])
		for Index in FaceIndices:
			Vertices = Vertices.union(Faces[InitialDim - 2][Index].Vertices)
		FaceSum = 0
		for i in xrange(len(Faces)):
			FaceSum += (-1)^(i+1)*len(Faces[i])
		PrintFaceLens()
		if len(Vertices) + FaceSum != (-1)^(InitialDim+1) + 1:
			print "Euler characteristic failed for polytope that lives in dimension", LocalDim
			print "len(Vertices) = ", len(Vertices)
			raw_input()
		else:
			print "Euler characteristic passed for polytope that lives in dimension", LocalDim
			return
	return FaceIndices

#-------------------------------------------------------------------------------
def MakeFace(FacePts, Pts, LongPointToShortPointMap, ShortPointToLongPointMap, Barycenter):
	Dimension = CheckNormalFormDim(GetHNF(FacePts))

	#This is to check if this set of points has already been found
	FacePtsInOriginalDim = []
	PtIndexSet = set()
	for Pt in FacePts:
		LongPt = ShortPointToLongPointMap[tuple(Pt)]
		FacePtsInOriginalDim.append(LongPt)
		PtIndexSet.add(PointToIndexMap[tuple(LongPt)])
	for i in xrange(len(Faces[Dimension-1])):
		if len((Faces[Dimension-1][i].Vertices).difference(PtIndexSet)) == 0:
			return i

	NewFace = Face()
	NewFaceIndex = len(Faces[Dimension-1])
	NewFace.Dimension = Dimension

	if Dimension == 2:
		FacePts, LocalShortPointToLongPointMap, LocalLongPointToShortPointMap, LocalDim = WrapMaps(FacePts)
		Hull = []
		for Pt in ConvexHull2d(FacePts):
			Hull.append(ShortPointToLongPointMap[tuple(LocalShortPointToLongPointMap[tuple(Pt)])])
		for Pt in Hull:
			NewFace.Vertices.add(PointToIndexMap[tuple(Pt)])

		for i in xrange(len(Hull)):
			EdgeFound = False
			Pt1 = PointToIndexMap[tuple(Hull[i])]
			Pt2 = PointToIndexMap[tuple(Hull[(i+1)%len(Hull)])]
			Edge = set([Pt1, Pt2])
			for i in xrange(len(Faces[0])):
				if Faces[0][i].Vertices == Edge:
					NewFace.Children.add(i)
					EdgeFound = True
					break
			if EdgeFound == False:
				EdgeFace = Face()
				EdgeFace.Dimension = 1
				EdgeFace.Vertices = Edge
				Faces[0].append(EdgeFace)
				NewFace.Children.add(len(Faces[0])-1)
				EdgeFace.Neighbors = set([])
				for i in xrange(len(Faces[0])):
					PossibleNeighbor = Faces[0][i]
					if len(PossibleNeighbor.Vertices.intersection(EdgeFace.Vertices)) == 1:
						PossibleNeighbor.Neighbors.add((len(Faces[0])-1,list(PossibleNeighbor.Vertices.intersection(EdgeFace.Vertices))[0]))
						EdgeFace.Neighbors.add((i,list(PossibleNeighbor.Vertices.intersection(EdgeFace.Vertices))[0]))
	elif Dimension > 2:
		NewFace.Children = GiftWrap(FacePtsInOriginalDim)
		for ChildFaceIndex in NewFace.Children:
			NewFace.Vertices = NewFace.Vertices.union(Faces[Dimension-2][ChildFaceIndex].Vertices)
	else:
		print "Internal error, dimension ", Dimension, " not expected in MakeFace."
		raw_input()

	FullDimVertices = []
	for Vertex in NewFace.Vertices:
		FullDimVertices.append(IndexToPointMap[Vertex])

	Normal = GetNormalFromHNF(GetHNF(FullDimVertices))
	NewFace.InnerNormal = MakeNormalPointInDirectionOfPt(Normal, Barycenter, FullDimVertices[0])

	for i in xrange(len(Faces[Dimension - 1])):
		PossibleNeighbor = Faces[Dimension - 1][i]
		Intersection = list(PossibleNeighbor.Children.intersection(NewFace.Children))
		if len(Intersection) == 1:
			PossibleNeighbor.Neighbors.add((NewFaceIndex,Intersection[0]))
			NewFace.Neighbors.add((i,Intersection[0]))

	Faces[Dimension-1].append(NewFace)
	IndexOfFace = NewFaceIndex
	return IndexOfFace


def MakeRandomPointSet(Dim,Num):
	Pts = []
	for i in xrange(Num):
		Pts.append([])
		for j in xrange(Dim):
			Pts[i].append(Integer(randint(-100,100)))
	return Pts

Tests = []
for i in xrange(1):
	Tests.append((MakeRandomPointSet(3,3), 'Random'))

"""
Cyclic = CreateCyclicLists(5)
Tests = []
for i in xrange(len(Cyclic)):
	Tests.append((Cyclic[i], 'Cyclic'))
"""

InitialTime = time()
for k in range(20):
	print ""
for Test in Tests:
	global Faces
	Faces = "Start"
	global PointToIndexMap
	PointToIndexMap = "Start"
	Pts = Test[0]
	print Test[1]
	print Pts
	StartTime = time()
	GiftWrap(Pts)
	print "Done with " + Test[1]
	print "Time = ", time() - StartTime
	print ""
	print ""
print "Total time = ", time() - InitialTime

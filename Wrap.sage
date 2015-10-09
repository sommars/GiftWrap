global Faces
Faces = "Start"
global PointToIndexMap
PointToIndexMap = "Start"
global InitialDim
InitialDim = "Sentinel"
global AllPts
load("GiftUtil.sage")
# The code in this unit is the code that will have to change for going between 3 and n dimensions.

#If all of the points have a column that's identical, this will break. I should immediately make a map in that case..?
def GiftWrap(Pts):
	global Faces
	global InitialDim
	global PointToIndexMap
	global IndexToPointMap
	global Barycenter
	global AllPts
	if not PtsAreValid(Pts):
		print "The input set of points is not valid."
		raw_input()
		return
	Pts = RemoveDups(Pts)
	OriginalBarycenter = FindBarycenter(Pts)
	if PointToIndexMap == "Start":
		AllPts = Pts
	Pts, ShortPointToLongPointMap, LongPointToShortPointMap, LocalDim = WrapMaps(AllPts, Pts)
	if PointToIndexMap == "Start":
		PointToIndexMap, IndexToPointMap = MakeIndexMaps(Pts)
		InitialDim = LocalDim
		AllPts = Pts
	if Faces == "Start":
		Faces = []
		for Dim in xrange(1,len(Pts[0])):
			Faces.append([])
	Barycenter = FindBarycenter(Pts)
	FirstPts = FindInitialFacet(Pts)
	IndexOfFirstFace = MakeFace(FirstPts, Pts, LongPointToShortPointMap, ShortPointToLongPointMap)
	FaceIndices = [IndexOfFirstFace]
	Counter = 0
	EdgesToRotateOver = set([])
	for ChildEdge in Faces[LocalDim-2][IndexOfFirstFace].Children:
		EdgesToRotateOver.add((IndexOfFirstFace, ChildEdge))
	while len(EdgesToRotateOver) > 0:
		if LocalDim == 4:
			print EdgesToRotateOver
		Edge = EdgesToRotateOver.pop()
		FaceWithEdge = Faces[LocalDim-2][Edge[0]]
		ChildFace = Faces[LocalDim-3][Edge[1]]


		Counter += 1
		if Counter == 200:
			print "Internal error, computation time exceeded"
			print len(Faces[0])
			print len(Faces[1])
			#PrintFaces()
			raw_input()

		VerticesToUse = []
		for Vertex in FaceWithEdge.Vertices:
			VerticesToUse.append(LongPointToShortPointMap[tuple(IndexToPointMap[Vertex])])


		EdgePts = []
		for Vertex in ChildFace.Vertices:
			EdgePts.append(LongPointToShortPointMap[tuple(IndexToPointMap[Vertex])])
		Normal = GetNormalFromHNF(GetHNF(VerticesToUse))
		# Need to make the normal an inner normal
		if not NormalPointsTowardsPt(Normal, Barycenter, VerticesToUse[0]):
			Normal = [-Normal[i] for i in xrange(len(Normal))]
		InnerNormal = Normal
		FacePts = EdgePts + FindNewFacetPtsFromEdge(Pts, EdgePts, InnerNormal, VerticesToUse)[0]
		#Call MakeFacet on the points I've found and do some bookkeeping
		NewFaceIndex = MakeFace(FacePts, Pts, LongPointToShortPointMap, ShortPointToLongPointMap)
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
		#PrintFaces()
		#for i in xrange(len(Faces[1])):
		#	print Faces[1][i].Vertices
		raw_input()
	#else:
		#print "Euler characteristic passed for polytope that lives in dimension", LocalDim
	if InitialDim == LocalDim:
		print "Euler characteristic passed for polytope that lives in dimension", LocalDim
		return 
	return FaceIndices

#-------------------------------------------------------------------------------
def Make2dFace(FacePts, Pts, LongPointToShortPointMap, ShortPointToLongPointMap):
	UCT = GetUCTAndNormal(GetNormalFromHNF(GetHNF(FacePts)))[0]
	ShiftedFacePts = TransformPts(FacePts, UCT)
	# I don't really need to remove points like this. I can just wait til the end and knock them off.
	Hull, PointsToRemove = ConvexHull2d(ShiftedFacePts)
	Hull = TransformPts(Hull, matrix(UCT^-1, ZZ))

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
	print IndexToPointMap
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
def MakeFace(FacePts, Pts, LongPointToShortPointMap, ShortPointToLongPointMap):
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
		return Make2dFace(FacePtsInCorrectDim, Pts, LongPointToShortPointMap, ShortPointToLongPointMap)
	elif Dimension > 2:
		NewFace = Face()
		NewFaceIndex = len(Faces[Dimension-1])
		NewFace.Children = set(GiftWrap(FacePtsInCorrectDim))
		NewFace.Vertices = set([])
		for ChildFaceIndex in NewFace.Children:
			ChildFace = Faces[Dimension-2][ChildFaceIndex]
			NewFace.Vertices = NewFace.Vertices.union(ChildFace.Vertices)
			ChildFace.Parents.append(NewFaceIndex)

		# the below lines should probably be refactored somehow since they show up multiple places
		# THIS IS GARBAGE. Evidence of bad design
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
		print "Internal error, dimension ", dimension, " not expected in MakeFace."
		raw_input()
	return IndexOfFace


for k in range(20):
	print ""

def MakeRandomPointSet(Dim):
	Pts = []
	for i in xrange(randint(20,20)):
		Pts.append([])
		for j in xrange(Dim):
			Pts[i].append(Integer(randint(0,100)))
	return Pts


#Tetrahedron = [[46, 64, 38], [68, 70, 19], [18, 97, 28], [20, 73, 66], [16, 81, 27], [37, 1, 16], [98, 25, 33], [29, 45, 45], [100, 19, 5], [83, 22, 68], [50, 95, 92], [62, 66, 33], [47, 48, 93], [82, 46, 41], [10, 13, 26], [53, 49, 51], [53, 74, 75]]
#Tests = [(Tetrahedron,'Tetrahedron')]
#GOTWEIRDBEHAVIORCHECKWHENWORKS[[22, 66, 45, 76], [45, 38, 6, 99], [59, 23, 22, 8], [8, 89, 10, 90], [67, 4, 4, 17], [73, 14, 50, 16], [27, 31, 18, 93], [48, 97, 12, 44], [80, 55, 56, 19], [4, 2, 34, 68], [65, 32, 7, 75], [61, 32, 2, 13], [48, 64, 8, 60], [61, 100, 91, 60], [41, 23, 85, 23], [80, 55, 48, 72], [31, 77, 2, 12], [60, 93, 16, 64], [53, 53, 50, 23], [55, 81, 9, 88]]

ZZZZ = [[47, 13, 48, 92], [8, 25, 69, 8], [85, 25, 75, 5], [86, 67, 96, 64], [95, 58, 89, 73], [17, 0, 35, 31], [100, 80, 81, 33], [63, 90, 11, 35], [71, 49, 85, 77], [60, 99, 25, 38], [66, 93, 44, 18], [84, 57, 58, 70], [19, 20, 76, 57], [1, 8, 87, 61], [82, 62, 53, 57], [82, 31, 8, 80], [51, 66, 0, 63], [46, 54, 37, 99], [35, 63, 65, 53], [9, 6, 60, 38]]
#ZZZZ = [[80, 51, 72], [100, 28, 36], [34, 49, 56], [26, 91, 86], [21, 55, 34], [41, 99, 63], [70, 34, 81], [48, 11, 78], [24, 94, 8], [52, 8, 81], [46, 99, 92], [99, 28, 90], [83, 2, 49], [27, 7, 87], [4, 20, 34], [99, 92, 76], [97, 4, 42], [16, 97, 53], [76, 70, 46], [45, 96, 45]]


#ZZZZ = [[93, 17, 5, 51], [34, 20, 11, 34], [100, 17, 86, 90], [37, 70, 86, 15], [49, 3, 28, 52], [91, 56, 15, 31], [59, 27, 32, 58], [53, 6, 84, 49], [5, 72, 60, 90], [47, 24, 46, 56], [86, 39, 26, 79], [98, 85, 23, 90], [56, 34, 12, 100], [89, 15, 7, 85], [40, 85, 39, 52], [69, 48, 97, 19], [60, 13, 95, 13], [29, 60, 92, 83], [73, 58, 18, 42], [82, 51, 23, 14]]

ZZZZ = [[47, 13, 48, 92], [8, 25, 69, 8], [85, 25, 75, 5], [86, 67, 96, 64], [95, 58, 89, 73], [17, 0, 35, 31]]#, [100, 80, 81, 33]]
#ZZZZ = [[80, 51, 72], [100, 28, 36], [34, 49, 56], [26, 91, 86], [21, 55, 34]]

Tests = []
for i in xrange(10):
	Tests.append((MakeRandomPointSet(3), 'Random'))

Tests = [(ZZZZ, 'ZZZZ')]

for Test in Tests:
	Pts = Test[0]
	print Test[1]
	print Pts
	GiftWrap(Pts)
	print "Done with " + Test[1]
	print ""
	print ""
	global Faces
	Faces = "Start"
	global PointToIndexMap
	PointToIndexMap = "Start"
	global InitialDim
	InitialDim = "Sentinel"

load("GiftUtil.sage")
# The code in this unit is the code that will have to change for going between 3 and n dimensions.
def GiftWrap(Pts):
	Facets = []
	if not PtsAreValid(Pts):
		return "The input set of points is not valid."
	Pts = RemoveDups(Pts)
	OriginalBarycenter = FindBarycenter(Pts)
	global InitialDim
	Pts, PointMap, InitialDim = MakePointMap(Pts)
	global Barycenter
	Barycenter = FindBarycenter(Pts)
	if InitialDim != 3:
		print "This set of point is not full dimensional"
		return
	FirstPts = FindInitialFacet(Pts, Barycenter)
	InitialFacet, PtsToRemove, EdgeList = MakeFacet(FirstPts, Pts)
	Facets.append(InitialFacet)
	Pts = RemovePts(Pts, PtsToRemove)
	# Going to spin through all of the edges in edgelist.
	# EdgeList should really be done via sets.
	Counter = 0
	while len(EdgeList) > 0:
		Edge = EdgeList[0]
		#Find Facet that has this Edge
		
		for Facet in Facets:
			Vertices = Facet.Vertices
			if (Edge[0] in Vertices) and (Edge[1] in Vertices):
				break

		Counter += 1
		if Counter == 200:
			return
			raw_input()
		FacetPts = Edge + FindNewFacetPtsFromEdge(Pts, Edge, Facet.InnerNormal, Facet.Vertices)

		#Call MakeFacet on the points I've found and do some bookkeeping
		Facet, PtsToRemove, Edges = MakeFacet(FacetPts, Pts)

		Pts = RemovePts(Pts, PtsToRemove)

		for Edge in Edges:
			# This is hacky. Really should just be careful and do some sorting
			ReversedEdge = [Edge[1], Edge[0]]
			if Edge in EdgeList:
				EdgeList.remove(Edge)
			elif ReversedEdge in EdgeList:
				EdgeList.remove(ReversedEdge)
			else:
				EdgeList.append(Edge)
		Facets.append(Facet)

	# If there are any points sitting in the interior of our polytope then they
	# have not yet been removed. We do so now. This is lazily programmed, and will
	# be fixed in the future.
	InteriorPts = list(Pts)
	for Facet in Facets:
		for Vertex in Facet.Vertices:
			if Vertex in InteriorPts:
				InteriorPts.remove(Vertex)
	Pts = RemovePts(Pts, InteriorPts)
	Facets = PutFacetsInCorrectDimension(Facets, PointMap, OriginalBarycenter)
	# Euler Characteristic check
	EdgeCount = 0
	for Facet in Facets:
		EdgeCount += len(Facet.Vertices)
	# We don't want to double count edges
	EdgeCount = EdgeCount/2
	if len(Pts) - EdgeCount + len(Facets) != 2:
		print "Euler Characteristic failed!"
		raw_input()
	# Note that in 4d it would be if len(Pts) - EdgeCount + len(3faces) - len(4faces) != 0:
	# In general, it's F_0 - F_1 + F_2 .... F_N != 1 + (-1)^(N+1)
	return Facets

#-------------------------------------------------------------------------------
def MakeFacet(FacetPts, Pts):
	# We have all of the points in their general position. We need to bring them
	# into vertical position so that the convexhull algorithm can work.
	UCT = GetUCTAndNormal(GetNormalFromHNF(GetHNF(FacetPts)))[0]
	ShiftedFacetPts = TransformPts(FacetPts, UCT)
	Hull, PointsToRemove = ConvexHull2d(ShiftedFacetPts)
	
	# We want to make sure the hyperplane supporting the facet is to one side
	# of all of the points.
	CheckAllPtsLieOnOthersideOfFacetHyperplane(copy(Pts), copy(ShiftedFacetPts), UCT)

	Hull = TransformPts(Hull, matrix(UCT^-1, ZZ))
	
	#Hate doing this. Not sure why it's necessary, but this is a temporary solution
	for i in xrange(len(Hull)):
		for j in xrange(len(Hull[i])):
			Hull[i][j] = Integer(Hull[i][j])

	PointsToRemove = TransformPts(PointsToRemove, matrix(UCT^-1, ZZ))
	Facet = Face()

	Normal = GetNormalFromHNF(GetHNF(Hull))
	# Need to make the normal an inner normal
	if not NormalPointsTowardsPt(Normal, Barycenter, Hull[0]):
		Normal = [-Normal[i] for i in xrange(len(Normal))]
	Facet.InnerNormal = Normal
	Edges = [[Hull[len(Hull) - 1], Hull[0]]]
	for i in xrange(len(Hull) - 1):
		Edges.append([Hull[i],Hull[i+1]])
	Facet.Vertices = Hull
	Facet.Children = Hull
	Dim = CheckNormalFormDim(GetHNF(Hull))
	if Dim != 2:
		print "Internal error, 2d convex hull is in the wrong dimension."
		raw_input()
	Facet.Dimension = Dim
	Facet.Neighbors = []
	return Facet, PointsToRemove, Edges
	
	
for k in range(20):
	print ""
TestCube = [[0,0,0],[1,0,0],[0,1,0],[1,1,0],[0,0,1],[1,0,1],[0,1,1],[1,1,1]]
BigCube = []
for i in [0,1,2]:
	for j in [0,1,2]:
		for k in [0,1,2]:
			BigCube.append([i,j,k])
TetrahedronExtraPts = [[1,1,1], [1,-1,-1], [-1,1,-1], [-1,-1,1],[0,0,0],[1,0,0]]
Tetrahedron = [[1,1,1], [1,-1,-1], [-1,1,-1], [-1,-1,1]]
TippedOverHouse = [[0,0,0],[2,0,0],[0,2,0],[2,2,0],[0,0,2],[2,0,2],[0,2,2],[2,2,2],[3,1,1]]
FiveDCube = [[1,0,0,0,0],[1,1,0,0,1],[1,0,1,0,0],[1,1,1,0,1],[1,0,0,1,0],[1,1,0,1,1],[1,0,1,1,0],[1,1,1,1,1]]

Tests = [(TestCube,'Cube'), (Tetrahedron,'Tetrahedron'),(TippedOverHouse,'TippedOverHouse'),(TetrahedronExtraPts,'TetrahedronExtraPts'),(BigCube,'BigCube')]

def MakeRandomPointSet():
	Pts = []
	for i in xrange(randint(20,20)):
		Pts.append([])
		for j in xrange(3):
			Pts[i].append(Integer(randint(0,100)))
	return Pts

Tests = []
for i in xrange(100):
	Tests.append((MakeRandomPointSet(), 'Random'))

#Tests = [(TestCube,'Cube'), (Tetrahedron,'Tetrahedron'),(TippedOverHouse,'TippedOverHouse'),(TetrahedronExtraPts,'TetrahedronExtraPts'),(BigCube,'BigCube')]
#Tests = [(FiveDCube,'3d Cube in 5d'), (TestCube,'Cube')]
#Tests = [(TestCube,'Cube')]
#Tests = [(TippedOverHouse,'TippedOverHouse')]#,(Tetrahedron,'Tetrahedron')]
#Tetrahedron = [[46, 64, 38], [68, 70, 19], [18, 97, 28], [20, 73, 66], [16, 81, 27], [37, 1, 16], [98, 25, 33], [29, 45, 45], [100, 19, 5], [83, 22, 68], [50, 95, 92], [62, 66, 33], [47, 48, 93], [82, 46, 41], [10, 13, 26], [53, 49, 51], [53, 74, 75]]
#Tests = [(Tetrahedron,'Tetrahedron')]

#ZZZZ = [[33, 36, 82], [7, 49, 16], [70, 95, 85], [41, 57, 43], [94, 99, 1], [91, 14, 69], [68, 6, 13], [33, 58, 53], [35, 24, 75], [75, 64, 79], [3, 80, 87]]
#ZZZZ = [[13, 18, 70], [68, 89, 61], [51, 1, 37], [37, 9, 88], [57, 56, 62], [77, 12, 58], [57, 4, 79], [39, 25, 33], [58, 57, 33], [59, 12, 94], [56, 35, 91], [90, 14, 38], [20, 19, 96], [87, 97, 15], [59, 49, 12], [14, 72, 4], [39, 93, 48]]



#ZZZZ = [[73, 50, 48], [7, 10, 66], [38, 62, 3], [38, 99, 3]]
#ZZZZ = [[59, 36, 7], [14, 95, 33], [70, 31, 33], [68, 6, 20]]
#ZZZZ = [[36, 31, 60, 94], [67, 29, 42, 63], [69, 94, 18, 61], [85, 30, 71, 100], [95, 2, 41, 99], [3, 28, 40, 79], [12, 83, 24, 22], [20, 53, 43, 95], [25, 81, 3, 84], [13, 85, 56, 61], [34, 76, 32, 55], [65, 22, 64, 52], [19, 39, 85, 72], [88, 71, 100, 95], [88, 85, 76, 59], [40, 41, 22, 87], [38, 2, 8, 91], [71, 43, 9, 50], [30, 8, 75, 62], [9, 68, 29, 14]]
#ZZZZ = [[65, 44, 65, 49], [38, 92, 47, 13], [33, 44, 95, 84], [56, 75, 7, 75], [98, 72, 95, 97], [82, 90, 4, 22], [79, 60, 86, 52], [52, 95, 33, 47], [17, 7, 98, 71], [18, 1, 70, 57], [37, 85, 52, 69], [95, 80, 65, 2], [93, 49, 77, 80], [66, 86, 29, 91], [75, 64, 12, 97], [100, 0, 15, 9], [99, 74, 98, 90], [42, 86, 36, 94], [55, 8, 21, 57], [73, 16, 34, 43]]
#ZZZZ = [[35, 82, 69, 28], [68, 52, 41, 90], [97, 66, 25, 58], [43, 63, 69, 5], [28, 14, 17, 1], [62, 48, 67, 91], [46, 73, 51, 25], [84, 55, 6, 81], [84, 43, 2, 58], [86, 62, 49, 21], [17, 76, 74, 16], [7, 20, 81, 57], [65, 9, 2, 4], [37, 44, 62, 5], [23, 65, 3, 78], [20, 65, 11, 67], [61, 16, 48, 13], [84, 64, 34, 74], [25, 77, 13, 26], [22, 61, 22, 17]]
#ZZZZ = [[91, 41, 53, 62], [78, 84, 59, 48], [37, 0, 17, 48], [88, 58, 34, 58], [22, 79, 20, 13], [80, 19, 33, 62], [94, 23, 9, 97], [85, 99, 79, 94], [93, 69, 57, 86], [17, 96, 21, 22], [87, 33, 5, 78], [68, 64, 33, 93], [45, 84, 53, 66], [98, 91, 99, 39], [95, 72, 65, 5], [14, 26, 83, 95], [11, 28, 84, 40], [56, 34, 50, 44], [78, 100, 35, 49], [49, 92, 56, 17]]
#ZZZZ = [[71, 31, 69, 68], [26, 76, 33, 71], [6, 67, 83, 30], [5, 79, 20, 3], [10, 10, 93, 88], [85, 31, 34, 5], [92, 28, 33, 91], [44, 97, 2, 19], [94, 0, 46, 39], [78, 88, 90, 37], [67, 98, 88, 61], [74, 75, 57, 39], [57, 62, 39, 7], [25, 41, 9, 93], [31, 21, 70, 19], [5, 69, 4, 68], [74, 11, 95, 2], [44, 19, 42, 22], [88, 12, 22, 30], [40, 95, 83, 64]]
#ZZZZ = [[53, 86, 37, 87], [7, 92, 13, 99], [96, 18, 61, 25], [62, 42, 81, 86], [11, 35, 89, 95], [65, 44, 75, 21], [55, 25, 71, 70], [100, 97, 99, 71], [31, 89, 95, 65], [1, 35, 24, 5], [51, 61, 6, 45], [100, 77, 81, 35], [90, 93, 58, 80], [82, 26, 49, 52], [99, 28, 95, 34], [68, 83, 8, 20], [40, 16, 64, 18], [15, 57, 27, 16], [19, 20, 43, 74], [58, 23, 32, 44]]
#Tests = [(ZZZZ, 'ZZZZ')]

for Test in Tests:
	Pts = Test[0]
	print Test[1]
	print Pts
	#FindInitialFacet(Pts)
	GiftWrap(Pts)
	print "Done with " + Test[1]
	print ""
	print ""

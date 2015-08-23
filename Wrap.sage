load("GiftUtil.sage")
# The code in this unit is the code that will have to change for going between 3 and n dimensions.
def GiftWrap(Pts):
	Facets = []
	if not PtsAreValid(Pts):
		return "The input is not valid."
	Pts = RemoveDups(Pts)
	global Barycenter
	Barycenter = FindBarycenter(Pts)
	# Currently, I'm going to blow up if my pointset isn't full dimensional
	if CheckNormalFormDim(GetHNF(Pts)) != 2:
		print "NOT FULL DIMENSITONAL"
		print Pts
		return
	InitialFacet, PtsToRemove, EdgeList = FindInitialFacetTwo(Pts)
	Facets.append(InitialFacet)
	Pts = RemovePts(Pts, PtsToRemove)
	# Going to spin through all of the edges in edgelist.
	# EdgeList should really be done via sets.
	Counter = 0
	DistinctEdges = list(EdgeList)
	while len(EdgeList) > 0:
		Edge = EdgeList[0]
		#Find Facet that has this Edge
		
		InnerNs = []
		for Facet in Facets:
			Facet.PrintProps()
			InnerNs.append(Facet.InnerNormal)
		print InnerNs

		for Facet in Facets:
			Vertices = Facet.Vertices
			if (Edge[0] in Vertices) and (Edge[1] in Vertices):
				break
		
		for k in xrange(5):
			print ""

		
		print Counter
		print EdgeList
		

	#	print "HI"
		Counter += 1
		if Counter == 10:
			raw_input()
		print "Edge", Edge
		print "FACET.Vertices", Facet.Vertices
		print "NEWPTS", FindNewFacetPtsTwo(Pts, Edge, Facet.InnerNormal, Facet.Vertices)
		FacetPts = Edge + FindNewFacetPtsTwo(Pts, Edge, Facet.InnerNormal, Facet.Vertices)

		#Call MakeFacet on the points I've found and do some bookkeeping
		Facet, PtsToRemove, Edges = MakeFacet(FacetPts)
		Pts = RemovePts(Pts, PtsToRemove)

		for Edge in Edges:
			Backward = [Edge[1], Edge[0]]
			if Edge not in DistinctEdges and Backward not in DistinctEdges:
				DistinctEdges.append(Edge)

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
	print "Finished! Count of faces = ", len(Facets)
	print ""
	print "....Computing Euler Characteristic...."
	print "Number of Vertices = ", len(Pts)
	EdgeCount = 0
	for Facet in Facets:
		#print Facet.Vertices
		EdgeCount += len(Facet.Vertices)
	EdgeCount = EdgeCount/2
	print "Number of Edges = ", EdgeCount
	print "Number of Facets", len(Facets)
	print len(Pts) - EdgeCount + len(Facets) == 2
	if len(Pts) - EdgeCount + len(Facets) != 2:
		print "Euler Characteristic failed!"
		raw_input()
	return Facets

#-------------------------------------------------------------------------------
def FindInitialFacetTwo(Pts):
	print "Start initial facet (version two)"
	def FindFirstPts(Pts):
		FirstPts = []
		FirstPtValue = Pts[0][0]
		for Pt in Pts:
			if Pt[0] == FirstPtValue:
				FirstPts.append(Pt)
			if Pt[0] > FirstPtValue:
				FirstPts = [Pt]
				FirstPtValue = Pt[0]
		return FirstPts
	FirstPts = FindFirstPts(Pts)

	if len(FirstPts) == 1:

		#Start with one point
		FirstPts = FirstPts + FindNewFacetPtsThree(Pts, FirstPts, [-1,0,0])
		
		#Now we have either an edge or a facet. Need to check which
		if CheckNormalFormDim(GetHNF(FirstPts)) == 1:
			print "Final FirstPtsList: ", FirstPts
			return MakeFacet(FirstPts)

		#if we get here, we have an edge. need to pick up one or more points
		#First, we need to have an inner normal
		InnerNormal = GetNormalFromHNF(GetHNF(FirstPts))
		if not NormalPointsTowardsPt(InnerNormal, Barycenter, FirstPts[0]):
			InnerNormal = [-InnerNormal[i] for i in xrange(len(InnerNormal))]

		KnownFacetPts = list(FirstPts)
		Vector = [FirstPts[0][i] - FirstPts[1][i] for i in xrange(len(FirstPts[0]))]
		Counter = 1
		while True:
			NewPoint = [FirstPts[0][i] + Counter*Vector[i] for i in xrange(len(FirstPts[0]))]
			if not NewPoint in KnownFacetPts:
				KnownFacetPts.append(NewPoint)
				break
			Counter += 1

		FirstPts = FirstPts + FindNewFacetPtsTwo(Pts, FirstPts, InnerNormal, KnownFacetPts)
		print "Final FirstPtsList: ", FirstPts
		return MakeFacet(FirstPts)		
	else:
		HNF = GetHNF(FirstPts)
		Dim = CheckNormalFormDim(HNF)
		if Dim == 1:
			print "Final FirstPtsList: ", FirstPts
			return MakeFacet(FirstPts)

#- However, if I have two points, this means that I have found an edge. I can do essentially the same thing I do in
#the general situation. I have the same normal to the hyperplane [-1,0,0]. I must use a specific normal to the pointset though, which I can #find by calling GetNormalFromHNF(GetHNFFromVectorAndPts([-1,0,0], Edge))
		elif Dim == 0:
			Normal = [-1,0,0]
			#Hacky
			KnownFacetPts = list(FirstPts)
			ExtraPt = [0 for i in xrange(len(FirstPts[0]))]
			ExtraPt[0] = FirstPts[0][0]
			for Pt in KnownFacetPts:
				for i in xrange(1,len(Pt)):
					ExtraPt[i] += abs(Pt[i])
			KnownFacetPts.append(ExtraPt)
			#EndHacky. Should refactor FindNewFacetPts routine so this isn't necessary
			FirstPts = FirstPts + FindNewFacetPtsTwo(Pts, FirstPts, [-1,0,0], KnownFacetPts)
			return MakeFacet(FirstPts)
		else:
			print "UNEXPECTED DIMENSION: ", Dim
			raw_input()
			return
	print "Final FirstPtsList: ", FirstPts
	return MakeFacet(FirstPts)

#-------------------------------------------------------------------------------
def MakeFacet(Pts):
	# We have all of the points in their general position. We need to bring them
	# into vertical position so that the convexhull algorithm can work.
	UCT = GetUCTAndNormal(GetNormalFromHNF(GetHNF(Pts)))[0]
	NewPts = TransformPts(Pts, UCT)
	Hull, PointsToRemove = ConvexHull2d(NewPts)

	Hull = TransformPts(Hull, matrix(UCT^-1, ZZ))
	#Hate doing this
	for i in xrange(len(Hull)):
		for j in xrange(len(Hull[i])):
			Hull[i][j] = Integer(Hull[i][j])
	PointsToRemove = TransformPts(PointsToRemove, matrix(UCT^-1, ZZ))
	Facet = Face()

	Normal = GetNormalFromHNF(GetHNF(Hull))
	# Need to check if the normal is actually an innernormal
	if not NormalPointsTowardsPt(Normal, Barycenter, Hull[0]):
		Normal = [-Normal[i] for i in xrange(len(Normal))]
	Facet.InnerNormal = Normal
	Hull.sort()
	Facet.Vertices = Hull
	Facet.Children = Hull
	Dim = CheckNormalFormDim(GetHNF(Hull))
	if Dim != 1:
		print "INCORRECT DIMENSION!!!"
		raw_input()
	Facet.Dimension = Dim
	Facet.Neighbors = []
	Edges = [[Hull[len(Hull) - 1], Hull[0]]]
	for i in xrange(len(Hull) - 1):
		Edges.append([Hull[i],Hull[i+1]])
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
Tests = [(TestCube,'Cube'), (Tetrahedron,'Tetrahedron'),(TippedOverHouse,'TippedOverHouse'),(TetrahedronExtraPts,'TetrahedronExtraPts'),(BigCube,'BigCube')]
#Tests = [(Tetrahedron,'Tetrahedron'),(TippedOverHouse,'TippedOverHouse'),(TetrahedronExtraPts,'TetrahedronExtraPts')]
#Tests = [(Tetrahedron,'Tetrahedron')]
#Tests = [(TetrahedronExtraPts,'TetrahedronExtraPts')]
#Tests = [(TippedOverHouse,'TippedOverHouse')]

def MakeRandomPointSet():
	Pts = []
	for i in xrange(randint(4,4)):
		Pts.append([])
		for j in xrange(3):
			Pts[i].append(Integer(randint(0,100)))
	return Pts
Tests = [(MakeRandomPointSet(), 'Random1')]

#Hangs indefinitely in the main loop
ZZZZ = [[81, 47, 67], [95, 63, 1], [85, 44, 84], [51, 62, 58], [73, 45, 65], [84, 38, 42], [48, 72, 91], [45, 89, 64], [9, 70, 78], [0, 59, 89], [39, 11, 70], [21, 69, 99], [28, 10, 18], [59, 40, 39], [20, 45, 19], [21, 73, 15], [55, 77, 93], [76, 54, 25], [96, 1, 68], [70, 77, 15], [42, 59, 50], [89, 33, 79], [11, 80, 74], [65, 37, 100], [5, 100, 73], [52, 54, 0], [0, 78, 68], [50, 67, 18], [36, 78, 68], [99, 37, 71], [65, 13, 31], [48, 51, 45], [91, 48, 43], [96, 13, 22], [30, 18, 38], [86, 29, 30], [23, 88, 94], [90, 87, 48], [64, 34, 44], [80, 93, 14], [4, 67, 82], [37, 55, 52], [21, 44, 36], [91, 95, 79]]

#Float division by zero issue
ZZZZ = [[72, 19, 88], [74, 40, 53], [45, 44, 58], [43, 92, 13], [10, 6, 16], [50, 45, 95], [33, 6, 59], [84, 74, 58], [88, 91, 20], [1, 100, 69], [35, 62, 44], [27, 57, 25], [76, 9, 78], [15, 50, 38], [31, 91, 19], [29, 65, 33], [77, 73, 88], [98, 3, 48], [30, 16, 1], [50, 93, 65], [41, 62, 33], [35, 77, 77], [57, 4, 22], [61, 83, 99], [81, 53, 35], [73, 37, 74], [99, 24, 92], [8, 64, 48], [73, 36, 5], [59, 60, 99], [94, 51, 12], [74, 50, 72], [0, 65, 20], [57, 92, 3], [68, 71, 77], [86, 5, 17], [33, 58, 98], [58, 54, 51], [49, 46, 47], [26, 54, 82], [14, 6, 25], [3, 66, 27], [57, 81, 15], [73, 0, 61], [19, 21, 54], [97, 30, 50], [8, 27, 49], [50, 29, 90], [3, 92, 77], [55, 64, 54], [70, 85, 17], [56, 40, 43], [18, 31, 67], [33, 21, 99], [75, 82, 50], [75, 33, 50], [34, 8, 19], [97, 75, 13], [71, 60, 95], [26, 26, 2], [46, 1, 83], [26, 75, 91], [72, 75, 55], [55, 98, 23], [74, 4, 12], [30, 18, 28], [82, 72, 48], [55, 90, 36], [12, 85, 41], [67, 55, 21], [49, 97, 99], [90, 73, 9], [95, 84, 18], [95, 7, 81], [29, 86, 35], [62, 48, 76], [94, 7, 98]]


#Hang in main loop
#ZZZZ = [[46, 48, 9], [99, 14, 6], [14, 91, 6], [77, 25, 81], [25, 72, 24], [56, 95, 46], [38, 2, 91], [56, 25, 79], [13, 96, 50], [43, 45, 22], [30, 89, 72]]

#ZZZZ = [[98, 15, 16], [11, 10, 99], [73, 84, 41], [93, 91, 69], [67, 85, 52], [25, 24, 49], [62, 77, 65], [56, 75, 33], [40, 33, 55], [74, 53, 53]]

ZZZZ = [[94, 28, 20], [92, 33, 93], [3, 81, 1], [73, 75, 15], [93, 90, 77]]


#Euler Characteristic failed
#ZZZZ = [[7, 99, 74], [10, 31, 38], [29, 7, 23], [61, 7, 72], [29, 91, 1], [19, 8, 27], [66, 9, 27], [51, 50, 2], [32, 96, 6], [17, 59, 0]]

Tests = [(ZZZZ, 'ZZZZ')]
Tests = []
for i in xrange(5000):
	Tests.append((MakeRandomPointSet(), 'Random1'))

Tests = [(ZZZZ, 'ZZZZ')]
for Test in Tests:
	Pts = Test[0]
	print Test[1]
	print Pts
#	print FindInitialFacetTwo(Pts)
	GiftWrap(Pts)
	print "Done with " + Test[1]
	print ""
	print ""

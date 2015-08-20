load("GiftUtil.sage")
# The code in this unit is the code that will have to change for going between 3 and n dimensions.
def GiftWrap(Pts):
	Facets = []
	if not PtsAreValid(Pts):
		return "The input is not valid."
	Pts = RemoveDups(Pts)
	global Barycenter
	Barycenter = FindBarycenter(Pts)
	
	InitialFacet, PtsToRemove, EdgeList = FindInitialFacet(Pts)
	Facets.append(InitialFacet)
	print "HIII"
	print PtsToRemove
	print Pts
	Pts = RemovePts(Pts, PtsToRemove)
	print Pts
	# Going to spin through all of the edges in edgelist.
	# EdgeList should really be done via sets.
	while len(EdgeList) > 0:
		Edge = EdgeList[0]

		#Find Facet that has this Edge
		for Facet in Facets:
			Vertices = Facet.Vertices
			if (Edge[0] in Vertices) and (Edge[1] in Vertices):
				break


		FacetPts = Edge + FindNewFacetPtsTwo(Pts, Edge, Facet.InnerNormal, Facet.Vertices)

		#Call MakeFacet on the points I've found and do some bookkeeping
		Facet, PtsToRemove, Edges = MakeFacet(FacetPts)
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
	print "FaceCount", len(Facets)
	print ""
	print "....Computing Euler Characteristic...."
	print "Number of Vertices = ", len(Pts)
	EdgeCount = 0
	for Facet in Facets:
		print Facet.Vertices
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
def FindInitialFacet(Pts):
	print "Start initial facet"
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
	print FirstPts
	Counter = 0
	while True:
		Counter += 1
		print Counter
		if len(FirstPts) != 1:	
			HNF = GetHNF(FirstPts)
			if CheckNormalFormDim(HNF) == 1:
				break
			Normal = GetNormalFromHNF(HNF)
		
		# If we're going through this the first time and we have a zero dimensional
		# set of points found so far, we want the normal to be the vector which we
		# minimized, i.e. [1,0,...,0], rather than what we get by taking the Hermite
		# Normal Form normal.
		if (Counter == 1):
			Normal = []
			Normal.append(1)
			for i in xrange(1,len(FirstPts[0])):
				Normal.append(0)
		UCT, NewNormal = GetUCTAndNormal(Normal)
		NewNormal = NewNormal[0]
		NewPts = TransformPts(Pts, UCT)
		NewFirstPts = TransformPts(FirstPts,UCT)

		# We check which way the normal should be pointed, positive or negative.
		# Note that I'm using outer normals. Another way of dealing with this would
		# be using the barycenter.
		if not NormalShouldBePositive(NewPts, NewFirstPts):
			NewNormal[0] = -1

		# Find a new point and add it to FirstPts. We add the points in the
		# original space since FirstPts is in the original space. Also, we transform
		# the normal to the new space, because otherwise its value is meaningless.
		FirstPts = FirstPts + TransformPts(GetAMaximalPt(NewPts, NewFirstPts, NewNormal), matrix(UCT^-1,ZZ))

		# Finally, we do another unimodular coordinate transform to find whatever 
		# other points lie on our facet that we haven't yet found.
		UCT = GetUCTAndNormal(GetNormalFromHNF(GetHNF(FirstPts)))[0]
		NewFirstPts = TransformPts(FirstPts, UCT)
		NewPts = TransformPts(Pts, UCT)
		VerticalPlaneCoord = NewFirstPts[0][0]
		for Pt in NewPts:
			if (Pt[0] == VerticalPlaneCoord) and (Pt not in NewFirstPts):
				NewFirstPts.append(Pt)
		FirstPts = TransformPts(NewFirstPts, matrix(UCT^-1, ZZ))

	print "Final FirstPtsList: ", FirstPts
	return MakeFacet(FirstPts)
"""
#-------------------------------------------------------------------------------
def FindInitialFacetTwo(Pts):
	print "Start initial facet"
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
	Counter = 0
	while True:
		Counter += 1
		if len(FirstPts) != 1:	
			HNF = GetHNF(FirstPts)
			if CheckNormalFormDim(HNF) == 1:
				break
			Normal = GetNormalFromHNF(HNF)
		
		# If we're going through this the first time and we have a zero dimensional
		# set of points found so far, we want the normal to be the vector which we
		# minimized, i.e. [1,0,...,0], rather than what we get by taking the Hermite
		# Normal Form normal.
		if (Counter == 1):
			Normal = []
			Normal.append(1)
			for i in xrange(1,len(FirstPts[0])):
				Normal.append(0)
		UCT, NewNormal = GetUCTAndNormal(Normal)
		NewNormal = NewNormal[0]
		NewPts = TransformPts(Pts, UCT)
		NewFirstPts = TransformPts(FirstPts,UCT)

		# We check which way the normal should be pointed, positive or negative.
		# Note that I'm using outer normals. Another way of dealing with this would
		# be using the barycenter.
		if not NormalShouldBePositive(NewPts, NewFirstPts):
			NewNormal[0] = -1

		# Find a new point and add it to FirstPts. We add the points in the
		# original space since FirstPts is in the original space. Also, we transform
		# the normal to the new space, because otherwise its value is meaningless.
		FirstPts = FirstPts + TransformPts(GetAMaximalPt(NewPts, NewFirstPts, NewNormal), matrix(UCT^-1,ZZ))
		
		# Finally, we do another unimodular coordinate transform to find whatever 
		# other points lie on our facet that we haven't yet found.
		UCT = GetUCTAndNormal(GetNormalFromHNF(GetHNF(FirstPts)))[0]
		NewFirstPts = TransformPts(FirstPts, UCT)
		NewPts = TransformPts(Pts, UCT)
		VerticalPlaneCoord = NewFirstPts[0][0]
		for Pt in NewPts:
			if (Pt[0] == VerticalPlaneCoord) and (Pt not in NewFirstPts):
				NewFirstPts.append(Pt)
		FirstPts = TransformPts(NewFirstPts, matrix(UCT^-1, ZZ))

	print "Final FirstPtsList: ", FirstPts
	return MakeFacet(FirstPts)
"""
#-------------------------------------------------------------------------------
def MakeFacet(Pts):
	# We have all of the points in their general position. We need to bring them
	# into vertical position so that the convexhull algorithm can work.
	UCT = GetUCTAndNormal(GetNormalFromHNF(GetHNF(Pts)))[0]
	NewPts = TransformPts(Pts, UCT)
	Hull, PointsToRemove = ConvexHull2d(NewPts)

	Hull = TransformPts(Hull, matrix(UCT^-1, ZZ))
	PointsToRemove = TransformPts(PointsToRemove, matrix(UCT^-1, ZZ))
	Facet = Face()

	Normal = GetNormalFromHNF(GetHNF(Hull))
	# Need to check if the normal is actually an innernormal
	if not NormalPointsTowardsPt(Normal, Barycenter, Hull[0]):
		Normal = [-Normal[i] for i in xrange(len(Normal))]
	Facet.InnerNormal = Normal

	Facet.Vertices = Hull
	Facet.Children = Hull
	Facet.Dimension = 2
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
	for i in xrange(randint(0,100)):
		Pts.append([])
		for j in xrange(3):
			Pts[i].append(Integer(randint(0,100)))
	return Pts
Tests = [(MakeRandomPointSet(), 'Random1')]
"""
for Test in Tests:
	Pts = Test[0]
	print Test[1]
	print Pts
#	print FindInitialFacetTwo(Pts)
	print GiftWrap(Pts)
	print "Done with " + Test[1]
	print ""
	print ""
"""
